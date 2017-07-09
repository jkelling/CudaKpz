/***************************************************************************
*   Copyright 2011 - 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifdef MODE_CPU
#include "kpzSchedulersCPU.h"
#else
#include "kpzSchedulersGPU.h"
#endif

#include <arghelper.h>
#include <dSFMT/dSFMT.h>
#include <benchmarking.h>
#include <timer.h>
#include <KMCsplash.h>
#include <kmcExceptCUDA.h>

#include <sstream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <limits>

#include "kpzFrontend.h"
using namespace Kpz;

void help()
{
	std::cout << "kpzScaling [options] mcs [x [y]]\n";
}

int main(int argc, char* argv[])
{
	timeAinit;
	registerTimeoutSignals();

	if(argc < 2)
	{
		help();
		return 1;
	}
	//get flags first
#ifdef MODE_OPENCL
	int kernelFlags = SchedulerService::KF_DEFAULT;
#endif
	std::string file;
	FrontendArgs args;
	{
		for(int a = args.parseSysArgs(argc,argv); a < argc; ++a)
		{
			if(!strcmp(argv[a], "-f") || !strcmp(argv[a], "--file"))
			{
				if(!getArg(file, ++a, argc, argv))
				{
					std::cerr << "Expected a filename for -f / --file.\n";
					return 1;
				}
			}
			else if(!strcmp(argv[a], "--help"))
			{
				help();
				version();
			}
#ifdef MODE_OPENCL
			else if(!strcmp(argv[a], "--innerDC"))
			{
				kernelFlags ^= SchedulerService::KF_INNER_DC;
			}
#endif
			else if((a = args.parseArg(argc,argv,a)) < 0)
				return 1;
		}
		args.init();
	}
#ifdef MODE_OPENCL
	SchedulerService* scheduler = new SchedulerService(x, y, kernelFlags);
#else
	SchedulerService* scheduler = args.schedulerConfig().mkScheduler(args.x(),args.x());
#endif
	// scheduler->supplySeed(23,42);

	std::cout << "mcsScale " << scheduler->mcsScale() << '\n';

	std::ofstream coutFile;
	if(file.length() && file != "--")
	{
		coutFile.open(file.c_str());
		if(!coutFile)
		{
			std::cerr << "Could not open file: '" << file << "' for writing.\n";
			return 1;
		}
		scheduler->setOut(coutFile);
	}
#ifdef MODE_CUDA
	if(!scheduler->size().validForGPU())
	{
		std::cerr << "System to small for current configuration.\n";
		return 1;
	}
#endif

	if(args.h5in() != 0)
	{
		auto data = scheduler->readH5(args.h5in());
		if(!data)
		{
			std::cerr << "[EE] Scheduler failed to read Hdf5 file.\n";
			exit(1);
		}
		try {
			unsigned int t;
			data->readGlobalAttribute("prevMeasureOverride", &t);
			args.pushPrevMeasureOverrideTime(t);
		}
		catch (splash::DCException) {}

		scheduler->closeH5(data);
		args.fastForwardEvents(scheduler->mcs());
	}
	else
	{
		timeAstart;
		scheduler->initHomogeneous();
		timeAstop("initSystem");
		timeAstart;
		scheduler->randomize();
		timeAstop("initRandomNumbers");
		if(!std::isnan(args.disorderFill()))
			scheduler->generateDisorder2(args.disorderFill());
	}

	void (SchedulerService::*mcsFkt) (int);
	if(args.mcsFct() == FrontendArgs::MCSFCT_NOP)
		mcsFkt = &SchedulerService::mcsNop;
	else if(!std::isnan(args.rateP()))
	{
		if(args.rateP() == 1. && args.rateQ() == 0.)
		{
			std::cerr << "[WW] Running dynamic pq implementation with p=1 and "
				"q=0, using the fixed implementation would be wore efficient in "
				"this case\n";
		}
		std::cout << "Rates: p = " << args.rateP() << " , q = " << args.rateQ() << '\n';
		if(!args.setPQ(scheduler))
		{
			std::cerr << "--disorder: setting parameters failed.\n";
			return 1;
		}
		if(std::isnan(args.disorderFill()))
			mcsFkt = &SchedulerService::mcsPQ;
		else
		{
			std::cout << "Disorder-rates: p1 = " << args.rateP(1) << " , q1 = " << args.rateQ(1) << '\n';
			mcsFkt = &SchedulerService::mcsDisorder2SyncRng;
		}
	}
	else
		mcsFkt = &SchedulerService::mcs;


	//scheduler->cout() << "0\t" << scheduler->roughness() << std::endl;
	// scheduler->roughnessAsync();
	scheduler->correlate();

	timeAstart;
	scheduler->pushSystem();
	timeAstop("pushSystem");

	unsigned int prevMeasureOverrideTime = scheduler->mcs();
	if(args.prevMeasureOverride() && args.overridePrevMeasurementTime(scheduler))
	{
		prevMeasureOverrideTime = args.prevMeasureOverrideTime();
	}

	while(!args.isFinished())
	{
		const unsigned int nextMcs = args.nextMcs(prevMeasureOverrideTime, scheduler);
		timeAstart;
		try {
			(scheduler->*mcsFkt)(nextMcs);
			timeAstopS( "mcs_" << nextMcs );
			scheduler->join();
			timeAstart;
			scheduler->popSystem();
			timeAstop("popSystem");
		}
		catch (KmcExcept::CUDAError e) {
			std::cout << "[EE] Caught CUDA error\n[EE] " << e.what() << "[EE] Going to timeout.\n";
			goto timeout;
		}
		if(Kmc::TimerSingleton::atEnd())
			goto timeout;
		scheduler->correlate();
		if(args.handleNextEvent(scheduler))
			prevMeasureOverrideTime = scheduler->mcs();
	}
	scheduler->join();

	if(false)
	{
timeout:
		scheduler->join();
		std::cout << "TIMEOUT " << scheduler->mcs() << '\n';
	}

	if(args.h5out() != 0)
	{
		auto data = scheduler->writeH5(args.h5out());
		data->writeGlobalAttribute(KmcSplash::ColTypeUInt32, "prevMeasureOverride", &prevMeasureOverrideTime);
		scheduler->closeH5(data);
	}
	if(args.fSave())
	{
		std::ostringstream ss;
		const char m_xyzPrefix[] = "KPZheightmap_";
			ss << m_xyzPrefix << scheduler->mcs() << ".xyz";
			scheduler->saveXYZ(ss.str().c_str());
	}

	delete scheduler;
	return 0;
}
