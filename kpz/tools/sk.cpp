#include <combineTables.h>
#include <arghelper.h>
#include <binnedVector.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <list>
#include <string>
#include <vector>
#include <cstring>
#include <algorithm>

static const size_t ALLOC = 10000;

struct Sample
{
	double val;
	// double sqval = NAN;
	unsigned int n;

		Sample(double val, unsigned int n)
			: val(val), n(n) {}
		Sample() : val(0), n(0) {}
};

struct StartTime
{
	unsigned int tmin;
	double k[4], avg, sumsq;
	double S, Q;
	unsigned long long count;
	BinnedVector<unsigned int> hist;
	BinnedVectorWindow histWindow {0,0};

	std::vector<Sample> samples;

		StartTime(unsigned int tmin)
			: tmin(tmin), S(0), Q(0), count(0), samples(0), avg(0.), sumsq(0.)
	{
		for(int i = 0; i < 4; ++i)
			k[i] = 0;
	}

	void initHist(unsigned int nBins, double yMin, double yMax) {
		histWindow = BinnedVectorWindow(yMin, yMax-yMin);
		hist.setBin((yMax-yMin)/nBins);
		histWindow.resize(hist, 0u);
		for(auto& a : hist)
			a = 0;
	}

	void printHist(std::ostream& o);

	static bool byTmin(const StartTime& a, const StartTime& b) {return a.tmin < b.tmin;}

	// calculate k* and count
	void eval(std::list<StartTime>::iterator beginOthers, std::list<StartTime>::iterator endOthers);

	void addSampleBatch_avg(const std::vector<Sample>& samples);
	void addSampleBatch_k(const std::vector<Sample>& samples);

	// calculate S and Q
	void calcSQ() {
		S = k[2] / pow(k[1],1.5);
		Q = k[3] / (k[1] * k[1]) - 3.0;
	}

	void calcSQAlt() {
		S = k[2] / (k[0] * k[0] * k[0]);
		Q = k[3] / (k[0] * k[0] * k[0] * k[0]);
	}

	double stdev() const {return sqrt((sumsq - avg*avg*count)/(count-1));}
};

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		std::cerr << "Usage: sk [options] [filenames] < [StartTimes with -t]\n";
		return 1;
	}
	std::list<int> idxY{1};
	// std::list<int> idxY2;
	int skipT = 0, skipOffset = 0;
	std::list<std::string> files;
	bool combined = false, enterTmin = false;
	std::list<StartTime> startTimes;
	bool intervals = false;
	unsigned int nBins = 200;
	std::string histPrefix("");

	for(int i = 1; i < argc; ++i)
	{
		if(!strcmp("-y", argv[i]))
		{
			std::string tmp;
			if(!getArg(tmp, ++i, argc, argv))
				return 1;
			if(!Kmc::parseIdxYList(idxY, tmp))
				return 1;
		}
		// else if(!strcmp("-y2", argv[i]))
		// {
		// 	std::string tmp;
		// 	if(!getArg(tmp, ++i, argc, argv))
		// 		return 1;
		// 	if(!Kmc::parseIdxYList(idxY2, tmp))
		// 		return 1;
		// }
		else if(!strcmp("-s", argv[i]))
		{
			if(!getArg(skipT, ++i, argc, argv))
				return 1;
			if(!getArg(skipOffset, ++i, argc, argv))
			{
				--i;
				skipOffset = 0;
			}
		}
		else if(!strcmp("-b", argv[i]))
		{
			if(!getArg(nBins, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("--histPrefix", argv[i]))
		{
			if(!getArg(histPrefix, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("-t", argv[i]))
		{
			enterTmin = true;
		}
		else if(!strcmp("--intervals", argv[i]))
		{ /* treat start times as intervalls, do not overlap runs. */
			intervals = true;
		}
		else if(!strcmp("-T", argv[i]))
		{
			unsigned int tmin;
			if(!getArg(tmin, ++i, argc, argv))
				return 1;
			startTimes.emplace_back(tmin);
		}
		else if(!strcmp("--combined", argv[i]))
		{
			combined = true;
		}
		else
		{
			files.push_back(argv[i]);
		}
	}

	const int nY = idxY.size();
	// const int nY2 = idxY2.size();
	// if(nY2 > 0)
	// {
	// 	if(nY2 != nY)
	// 	{
	// 		std::cerr << "[EE] #idxY2 must be equal to #idxY.\n";
	// 		return 1;
	// 	}
	// 	if(std::max(idxY) >= std::min(idxY2))
	// 	{
	// 		std::cerr << "[EE] All idxY2 must be greater than each idxY.\n";
	// 		return 1;
	// 	}
	// 	idxY.splice(idxY.end(), idxY2);
	// }

	{
		unsigned int tmin;
		if(enterTmin)
			while(std::cin >> tmin)
				startTimes.emplace_back(tmin);
		if(startTimes.empty())
		{
			tmin = 1000000;
			for(int a = 1; a < 10; ++a, tmin+=1000000)
				startTimes.emplace_back(tmin);
			for(int a = 1; a < 10; ++a, tmin+=10000000)
				startTimes.emplace_back(tmin);
		}
		startTimes.sort(StartTime::byTmin);
	}

	/* colX = 0, colY = idxY, colN = 4 (if combined) */
	double yMin = INFINITY, yMax = -INFINITY;
	unsigned int cnt = 0;
	for(auto file = files.begin(); file != files.end(); ++file)
	{
		Kmc::TableFileView table(0,idxY, file->c_str(), (combined)?4:-1);
		if(!table.open())
		{
			std::cerr << "[WW] Cannot read file: " << *file << '\n';
			continue;
		}

		auto currentStartTime = startTimes.begin();
		auto nextStartTime = currentStartTime;
		++nextStartTime;

		while(table.good() && table.x() < currentStartTime->tmin)
			table.nextLine();

		for(int nextTime = table.x()+skipOffset;table.good(); table.nextLine())
		{
			if(table.x() < nextTime)
				continue;
			nextTime = table.x() + skipT;
			if(nextStartTime != startTimes.end() && table.x() >= nextStartTime->tmin)
			{
				currentStartTime = nextStartTime;
				++nextStartTime;
			}
			if(currentStartTime->samples.capacity() == currentStartTime->samples.size())
			{
				currentStartTime->samples.reserve(currentStartTime->samples.size()+ALLOC);
			}

			for(int a = 0; a < nY; ++a)
			{
				if(a < yMin)
					yMin = a;
				if(a > yMax)
					yMax = a;
				// if(nY2 == 0)
				// 	currentStartTime->samples.emplace_back(table.y(a), table.n(), table.y(a+nY));
				// else
				currentStartTime->samples.emplace_back(a, table.n());
				++cnt;
			}
		}
	}
	if(cnt < 2)
	{
		std::cerr << "[EE] Got only " << cnt << " y-entries. Aborting.\n[MM] try adjusting startTimes.\n";
		return 2;
	}

	for(auto a = startTimes.begin(); a != startTimes.end(); ++a)
	{
		a->initHist(nBins, yMin, yMax);
		if(intervals)
			a->eval(a, a);
		else
		{
			auto b = a;
			a->eval(++b, startTimes.end());
		}
	}

	if(skipT != 0 || skipOffset != 0)
		std::cout << "#skip " << skipT << " +" << skipOffset << '\n';

	std::cout << "#tmin\tS\tQ\tcount\tavg\tstdev\tstderr\n";
	for(auto a = startTimes.begin(); a != startTimes.end(); ++a)
	{
		a->calcSQ();
		const auto stdev = a->stdev();
		std::cout << a->tmin << '\t' << a->S << '\t' << a->Q << '\t' << a->count
			<< '\t' << a->avg << '\t' << stdev << '\t' << stdev/sqrt((double)a->count) << '\n';

		if(!histPrefix.empty())
		{
			std::ostringstream os;
			os << histPrefix << '_' << a->tmin << ".dat";
			std::ofstream histOut(os.str().c_str());
			a->printHist(histOut);
		}
	}

	// std::cout << "#alternative rule\n";
	// for(auto a = startTimes.begin(); a != startTimes.end(); ++a)
	// {
	// 	a->calcSQAlt();
	// 	std::cout << a->tmin << '\t' << a->S << '\t' << a->Q << '\t' << a->count << '\n';
	// }

	return 0;
}

void StartTime::eval(std::list<StartTime>::iterator beginOthers, std::list<StartTime>::iterator endOthers)
{
	addSampleBatch_avg(samples);
	for(auto a = beginOthers; a != endOthers; ++a)
		addSampleBatch_avg(a->samples);

	avg /= count;

	addSampleBatch_k(samples);
	for(auto a = beginOthers; a != endOthers; ++a)
		addSampleBatch_k(a->samples);

	for(int i = 0; i < 4; ++i)
		k[i] /= count;
}

void StartTime::addSampleBatch_avg(const std::vector<Sample>& samples)
{
	for(auto a = samples.begin(); a != samples.end(); ++a)
	{
		count += a->n;
		const double add = a->n*a->val;
		avg += add;
		sumsq += add*a->val;

		histWindow(hist, a->val) += a->n;
	}
}

void StartTime::addSampleBatch_k(const std::vector<Sample>& samples)
{
	for(auto a = samples.begin(); a != samples.end(); ++a)
	{
		double d = a->val-avg;
		double acc = d;
		for(int i = 0; i < 4; ++i, acc *= d)
			k[i] += acc*a->n;
	}
}

void StartTime::printHist(std::ostream& o)
{
	for(int a = 0; a < hist.size(); ++a)
	{
		const double val = hist.at(a)/(double)count/hist.bin();
		o << histWindow.min()+a*hist.bin() << '\t' << val << '\t';
		o << (histWindow.min()+a*hist.bin())/avg << '\t' << val*avg << '\n';
	}
}
