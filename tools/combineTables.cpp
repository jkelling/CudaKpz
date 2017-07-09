/***************************************************************************
*   Copyright 2011 - 2016 Jeffrey Kelling <j.kelling@hzdr.de>
*                  Helmholtz-Zentrum Dresden-Rossendorf
*                  Institute of Ion Beam Physics and Materials Research
*
*	This file is part of CudaKpz.
*
*   CudaKPZ is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   CudaKPZ is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with CudaKPZ.  If not, see <http://www.gnu.org/licenses/>.
***************************************************************************/

#include <limits>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <memory>
#include <utility>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <glob.h>

#include <arghelper.h>
#include <combineTables.h>

using namespace Kmc;

class WhiteListCheckerAllowAll
{
	public:
		virtual bool atEnd() const {return false;}
		virtual bool check(const double x) {return true;}
		virtual void reset() {}

		virtual ~WhiteListCheckerAllowAll() {}
};

class WhiteListChecker : public WhiteListCheckerAllowAll
{
	std::vector<double>::const_iterator m_current;
	const std::vector<double>& m_list;

	public:

		WhiteListChecker(const std::vector<double>& list)
			: m_current(list.begin()), m_list(list)
		{}

		virtual bool atEnd() const {return m_list.end() == m_current;}
		virtual bool check(const double x) {
			while((!atEnd()) && x > *m_current)
			{
				WARN() << "[WW] skipping xWhitelist-entry: " << *m_current << '\n';
				++m_current;
			}
			if(atEnd())
				return false;
			if(x == *m_current)
			{
				++m_current;
				return true;
			}
			else return false;
		}
		virtual void reset() {m_current = m_list.begin();}
};

template<class Sample = Kmc::Sample>
std::list<Sample> processAvg(std::ostream& out, std::list<std::pair<TableFileView*, double> > files
	, unsigned int minSamples, unsigned int maxSamples, const std::string& metaKey, const std::vector<double>& xWhitelist
	, long long xmin = std::numeric_limits<long long>::min() , long long xmax = std::numeric_limits<long long>::max())
{
	std::list<Sample> samples;
	std::unique_ptr<WhiteListCheckerAllowAll> wlc;
	if(xWhitelist.empty())
		wlc.reset(new WhiteListCheckerAllowAll());
	else
		wlc.reset(new WhiteListChecker(xWhitelist));
	for(auto a = files.begin(); a != files.end(); ++a)
	{
		if(!a->first->open())
		{
			std::cerr << "[EE] Failed to open file " << a->first->name() << '\n';
			continue;
		}
		auto s = samples.begin();
		wlc->reset();

		auto addLine = [&]() -> void {
			const long long x = a->first->x();
			for(; s != samples.end() && s->x() < x; ++s);
			if(s != samples.end() && s->x() == x)
			{
				if(s->n() < maxSamples)
				{
					if(a->second != 1.)
						s->addWeighted(*a->first, a->second);
					else
						(*s) += *a->first;
				}
			}
			else
			{
				s = samples.insert(s, Sample(x, a->first));
				if(a->second != 1.)
					s->addWeighted(*a->first, a->second);
				else
					(*s) += *a->first;
			}
		};

		if(!metaKey.empty())
			while((a->first->meta() != metaKey || !wlc->check(a->first->x())) && a->first->good())
				a->first->nextLine();
		if(a->first->good())
		{
			while((a->first->x() <= xmin) && a->first->nextLine().good());
			if(a->first->x() > xmax)
				continue;
			addLine();
		}
		while(a->first->nextLine().good())
		{
			if(a->first->x() > xmax)
				break;
			if(a->first->x() < xmin)
				continue;
			if(!metaKey.empty())
			{
				while(a->first->meta() != metaKey && a->first->nextLine().good());
				if(!a->first->good())
					break;
			}
			if(a->first->x() <= s->x()) // skip double line
				continue;
			if(!wlc->check(a->first->x()))
			{
				if(wlc->atEnd())
					break; // nothing left in whitelist, discard rest
				else
					continue;
			}
			++s;
			addLine();
		}
		a->first->close();
	}
	out << "#x\tavg\tstdev\tsterr\tsamples\n";
	for(auto a = samples.begin(); a != samples.end(); ++a)
	{
		if(a->n() >= minSamples)
			out << (*a) << '\n';
	}

	return std::move(samples);
}

int main (int argc, char* argv[])
{
	int idxX = 0, idxN = -1, idxSigma = -1, idxMeta = -1;
	unsigned int minSamples = 0, maxSamples = std::numeric_limits<unsigned int>::max();
	std::list<int> idxY{1};
	std::list<std::pair<TableFileView*, double> > files;
	std::string metaKey;
	std::vector<double> xWhitelist;
	double weight = 1.;
	unsigned int hist = 0;
	double binWidth = NAN;
	long long xmax = std::numeric_limits<long long>::max();
	long long xmin = std::numeric_limits<long long>::min();
	const char* outFile = 0;
	
	for(int i = 1; i < argc; ++i)
	{
		if(!strcmp("-x", argv[i]))
		{
			if(!getArg(idxX, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("--minSamples", argv[i]))
		{
			if(!getArg(minSamples, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("--hist", argv[i]))
		{
			if(!getArg(hist, ++i, argc, argv))
			{
				hist = 200;
			}
		}
		else if(!strcmp("--bin", argv[i]))
		{
			hist = 1;
			if(!getArg(binWidth, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("--xmax", argv[i]))
		{
			if(!getArg(xmax, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("--xmin", argv[i]))
		{
			if(!getArg(xmin, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("--out", argv[i]))
		{
			if(!getArg(outFile, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("--maxSamples", argv[i]))
		{
			if(!getArg(maxSamples, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("-y", argv[i]))
		{
			std::string tmp;
			if(!getArg(tmp, ++i, argc, argv))
				return 1;
			if(!parseIdxYList(idxY, tmp))
				return 1;
		}
		else if(!strcmp("-n", argv[i]))
		{
			if(!getArg(idxN, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("-s", argv[i]))
		{
			if(!getArg(idxSigma, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("-p", argv[i]))
		{
			unsigned int tmp;
			if(!getArg(tmp, ++i, argc, argv))
				return 1;
			std::cout << std::setprecision(tmp);
		}
		else if(!strcmp("--select", argv[i]))
		{
			if(!metaKey.empty())
			{
				std::cerr << "[EE] Can only set one --select.\n";
				return 1;
			}
			if(!getArg(idxMeta, ++i, argc, argv))
				return 1;
			if(!getArg(metaKey, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("-w", argv[i]))
		{
			if(!getArg(weight, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("--combined", argv[i]))
		{
			idxX = 0;
			idxN = 4;
			idxSigma = 2;
			idxMeta = -1;
			idxY.clear();
			idxY.push_back(1);
			metaKey.clear();
		}
		else if(!strcmp("-v", argv[i]) || !strcmp("--verbose", argv[i]))
		{
			WARN.on();
		}
		else if(!strcmp("--xWhitelist", argv[i]))
		{
			const char* tmp;
			if(!getArg(tmp, ++i, argc, argv))
				return 1;
			std::vector<std::vector<double> > dummy;
			TableFileView xWhitelistFile(idxX, tmp, idxMeta);
			xWhitelistFile.read(xWhitelist, dummy);
			// for(auto i : xWhitelist)
			// 	std::cerr << i << '\n';
			// return 0;
		}
		else
		{
			glob_t names;
			if(glob(argv[i], GLOB_MARK, 0, &names))
				WARN() << "No files for pattern '" << argv[i] << "'\n";
			else
			{
				for(size_t a = 0; a < names.gl_pathc; ++a)
				{
					if(names.gl_pathv[a][strlen(names.gl_pathv[a])-1] == '/')
						std::cerr << "[EE] is a directory '" << names.gl_pathv[a] << "'\n";
					else
						files.push_back(
							std::make_pair(new TableFileView(idxX, idxY, names.gl_pathv[a], idxN, idxSigma, idxMeta), weight));
				}
			}
			globfree(&names);
		}
	}

	if(hist > 0)
	{
		auto samples = std::move(processAvg<CacheSample>(std::cout, files, minSamples, maxSamples, metaKey, xWhitelist, xmin, xmax));
		for(auto& s : samples)
		{
			if(!s.histGood())
				continue;
			std::ostringstream os;
			os << "hist_t_" << s.x() << ".dat";
			std::ofstream of(os.str().c_str());
			if(std::isnan(binWidth))
				s.printHist(of, hist);
			else
				s.printHist(of, s.nbinsFromWidth(binWidth));
		}
	}
	else
	{
		if(!outFile || !strcmp(outFile, "--"))
			processAvg(std::cout, files, minSamples, maxSamples, metaKey, xWhitelist);
		else
		{
			std::ofstream of(outFile);
			if(!of)
			{
				std::cerr << "[EE] failed to open outFile: " << outFile << '\n';
			}
			processAvg(of, files, minSamples, maxSamples, metaKey, xWhitelist);
		}
	}

	for(auto a = files.begin(); a != files.end(); ++a)
	{
		// std::cerr << *(*a);
		delete a->first;
	}

	return 0;
}
