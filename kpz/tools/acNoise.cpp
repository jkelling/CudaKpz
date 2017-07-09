#include <combineTables.h>
#include <arghelper.h>

#include <iostream>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstring>

struct Deviation
{
	double lamcz, s, c;
	
	Deviation(double lamcz, double s, double c=0)
		: lamcz(lamcz), s(s), c(c) {}

	double operator() (double x, double y) const {
		return y*(pow(x/s,lamcz)) - c;
	}
};

std::ostream& operator<<(std::ostream& o, const Deviation& d) {
	return o << "lamcz= " << d.lamcz << " ,s= " << d.s << " ,c= " << d.c;
}

inline double sq(double a) {
	return a*a;
}

int main(int argc, char* argv[])
{
	Deviation deviation(2.35, 100.);
	int idx = 0;
	int idy = 6;
	int minSamples = 10;
	unsigned int tmax = std::numeric_limits<unsigned int>::max();
	double tsMin = 2;
	std::ostream& out = std::cout;

	std::list<char*> files;
	for(int i = 1; i < argc; ++i)
	{
		if(!strcmp("-x", argv[i]))
		{
			if(!getArg(idx, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("-y", argv[i]))
		{
			if(!getArg(idx, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("-ms", argv[i]))
		{
			if(!getArg(minSamples, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("--fit", argv[i]))
		{
			if(!getArg(deviation.lamcz, ++i, argc, argv))
				return 1;
			if(!getArg(deviation.s, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("--tsMin", argv[i]))
		{
			if(!getArg(tsMin, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("--tmax", argv[i]))
		{
			if(!getArg(tmax, ++i, argc, argv))
				return 1;
		}
		else
		{
			files.push_back(argv[i]);
		}
	}

	std::map<unsigned int, std::vector<double> > data;

	{
		int filesLeft = files.size();
		for(auto a = files.begin(); a != files.end(); ++a, --filesLeft)
		{
			Kmc::TableFileView view(idx, idy, *a);
			if(!view.open())
			{
				std::cerr << "[EE] Failed to open file " << view.name() << '\n';
				continue;
			}

			for(;view.good(); view.nextLine())
			{
				const unsigned int x(view.x());
				if(x/deviation.s < tsMin)
					continue;
				if(x > tmax)
					break;
				auto t = data.find(x);
				if(t != data.end())
					t->second.push_back(view.y(0));
				else
				{
					auto &v = data[x];
					v.reserve(filesLeft);
					v.push_back(view.y(0));
				}
			}
		}
	}

	unsigned int nSamples = std::numeric_limits<unsigned int>::max();
	std::list<std::pair<unsigned int, std::vector<double> > > table;
	for(auto a = data.begin(); a!=data.end();)
	{
		if(a->second.size() >= minSamples)
		{
			nSamples = std::min((unsigned int)a->second.size(), nSamples);
			table.push_back(std::make_pair(a->first, std::vector<double>()));
			++a;
		}
		else
			a = data.erase(a);
	}
	if(nSamples == std::numeric_limits<unsigned int>::max())
	{
		std::cerr << "[EE] Did not find enough samples for any line to fulfill minSamples criterion (>= " << minSamples << ")\n";
		return 1;
	}
	out << "# truncating to " << nSamples << " samples\n";
	table.sort([](const std::pair<double, std::vector<double> >& a, const std::pair<double, std::vector<double> >& b){
			return a.first<b.first;});

	unsigned int maxRlSqrt = sqrt(nSamples);
	for(auto l = table.begin(); l != table.end(); ++l)
		l->second.resize(maxRlSqrt+1);

	// do runLength = nSamples to abtain deviation.c
	{
		double cavg = 0.;
		for(auto l = table.begin(); l != table.end(); ++l)
		{
			const auto& line = data[l->first];
			double avg = 0.;
			for(int a = 0; a < nSamples; ++a)
				avg += line[a];
			avg /= nSamples;
			cavg += l->second[0] = deviation(l->first, avg);
		}
		cavg /= table.size();
		deviation.c = cavg;
		for(auto l = table.begin(); l != table.end(); ++l)
			l->second[0] = fabs(l->second[0] - cavg);
	}
	out << "# " << deviation << '\n';

	for(auto l = table.begin(); l != table.end(); ++l)
	{
		const auto& line = data[l->first];
		for(int runLengthSqrt = 1; runLengthSqrt <= maxRlSqrt; ++runLengthSqrt)
		{
			const auto runLength = runLengthSqrt*runLengthSqrt;
			int nRuns = 0;
			l->second.at(runLengthSqrt) = 0.;
			for(int a = 0; a+runLength <= nSamples;)
			{
				double avg = 0.;
				for(int s = 0; s < runLength; ++s, ++a)
					avg += line[a];
				avg /= runLength;
				++nRuns;
				l->second[runLengthSqrt] += sq(deviation(l->first, avg));
			}
			l->second[runLengthSqrt] = sqrt(l->second[runLengthSqrt] / nRuns);
		}
	}

	// print
	out << "#t\t" << nSamples;
	for(int runLengthSqrt = 1; runLengthSqrt <= maxRlSqrt; ++runLengthSqrt)
		out << '\t' << runLengthSqrt*runLengthSqrt;
	out << '\n';

	for(auto l = table.begin(); l != table.end(); ++l)
	{
		out << l->first;
		for(auto a = l->second.cbegin(); a != l->second.cend(); ++a)
			out << '\t' << *a;
		out << '\n';
	}

	return 0;
}
