/***************************************************************************
*   Copyright 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#include "schedulerService.h"
#include "correlator.h"

#include <arghelper.h>
#include <benchmarking.h>
#include <KMCsplash.h>

#include <sstream>
#include <fstream>
#include <regex>
#include <cstring>

using namespace Kpz;

void help()
{
	std::cout << "kpzConsistencyCheck [--anaThreads num] [files]\n";
}

int main(int argc, char* argv[])
{
	// std::regex RX_filename("(.*)/(.*)_([01])_0_0.h5"); // subexressions are: dirname, name, mpiRank
	std::regex RX_filename("(.*/)(.*)_(.)_0_0\\.h5", std::regex::awk); // subexressions are: dirname, name, mpiRank
	for(int a = 1; a < argc; ++a)
	{
		if(!strcmp(argv[a], "--anaThreads"))
		{
			int n;
			if(!getArg(n, ++a, argc, argv))
			{
				std::cerr << "Expected number of threads for --anaThreads.\n";
				return 1;
			}
			std::cout << "[MM] running " << n << " threads for analysis (per process)\n";
			Kpz::Correlator::setNThreads(n);
		}
		else
		{
			std::string fname(argv[a]);
			std::smatch match;
			if(!std::regex_match(fname, match, RX_filename))
			{
				std::cerr << "[WW] Skipping file with non-matching name: " << fname << '\n';
				continue;
			}
			SchedulerService scheduler(4,4);
			std::string prefix;
			if( match[3] == "1" )
				prefix = "dis"; // assuming disordered part of autoResponse
			if(!scheduler.readH5(fname.c_str(), prefix))
			{
				std::cerr << "[EE] Failed to read file: " << fname << '\n';
				continue;
			}
			if( scheduler.consitencyCheck(std::cout) )
			{
				std::cout << "[MM] Good: " << fname << '\n';
			}
			else
				std::cout << "[MM] Bad: " << fname << '\n';
		}
	}
	return 0;
}
