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

#include "kpzFrontend.h"

#include "kpzConst.h"
#include "correlator.h"

#include <arghelper.h>
#include <timer.h>

#include <cstring>
#include <utility>
#include <algorithm>
#include <limits>
#include <iomanip>
#include <signal.h>

void Kpz::version()
{
	std::cout << "kpz* [vNAN]"
		"\nKMC_SWITCH_BENCHMARK "
#ifdef KMC_SWITCH_BENCHMARK
		"on"
#else
		"off"
#endif
		"\nKPZ_SWITCH_RNG "
#define xstr(s) str(s)
#define str(s) #s
		xstr(KPZ_SWITCH_RNG)
#undef xstr
#undef str
		"\nKPZ_ASYNC_MEASUREMENT_CPU "
#ifdef KPZ_ASYNC_MEASUREMENT_CPU
		"on"
#else
		"off"
#endif
#ifdef KPZ_FERMI
		"\n -- compiled for Fermi architecture"
#elif defined KPZ_K80
		"\n -- compiled for Kepler GK210 architecture"
#endif
		"\n";
}

void Kpz::registerTimeoutSignals()
{
	signal(SIGUSR1, sigINTtimeout);
	signal(SIGINT, sigINTtimeout);
}

void Kpz::sigINTtimeout (int sig)
{
	std::cerr << "Cought SIGINT: timeout.\n";
	Kmc::TimerSingleton::create((time_t)0);
}

Kpz::SchedulerConfig::SchedulerConfig()
	: m_schedulerId("default")
	  , m_db_lORest(-1)
	  , m_gpuYield(false), m_fAnalysisOnGPU(false)
	  , m_localLayoutPreset("default")
{
	m_ddBlocksize[0] = 0;
	m_ddBlocksize[1] = 0;
}

void Kpz::SchedulerConfig::setFlags(SchedulerService* scheduler) const
{
	scheduler->setAnalysisOnGPU(m_fAnalysisOnGPU);
}

int Kpz::SchedulerConfig::parseArgs(int argc, char* argv[], int a, int preferredSchedulerType)
{
	m_preferredSchedulerType = preferredSchedulerType;
	if(!getArg(m_schedulerId, a, argc, argv))
	{
		std::cerr << "Expected scheduler id.\n";
		return -1;
	}
	if(m_schedulerId == "+db")
	{
		m_schedulerId = "db";
		std::string tmp;
		if(getArg(tmp, ++a, argc, argv))
		{
			if(tmp == "-r")
			{
				// origin restriction, default: 0 (r1)
				if(!getArg(m_db_lORest, ++a, argc, argv))
				{
					std::cerr << "Expected integer arg for \\db -r.\n";
					return -1;
				}
			}
			else
				a--;
		}
	}

	if(m_schedulerId == "sca");
	else if(m_schedulerId == "default"); // dt or random sequential (cpu)
	else if(m_schedulerId == "dt");
	else if(m_schedulerId == "db4") // GPU only
	{
		if(m_preferredSchedulerType == SCHEDULER_TYPE_CPU)
			std::cerr << "[WW][SchedulerConfig] Scheduler " << m_schedulerId << " not available for CPU, will fall back to GPU.\n";
		m_preferredSchedulerType = SCHEDULER_TYPE_GPU;
	}
	else if(m_schedulerId == "db" || (m_schedulerId == "dtStatic")) // CPU only
	{
		if(m_preferredSchedulerType == SCHEDULER_TYPE_GPU)
			std::cerr << "[WW][SchedulerConfig] Scheduler " << m_schedulerId << " not available for GPU, will fall back to CPU.\n";
		m_preferredSchedulerType = SCHEDULER_TYPE_CPU;
	}
	else if(m_schedulerId == "scaBit" ); // SchedulerCPUCBBit or SchedulerSCABit
	else if(m_schedulerId == "cb")
		m_schedulerId = "sca";
	else
	{
		std::cerr << "[EE][SchedulerConfig] No such scheduler: " << m_schedulerId << '\n';
		return -1;
	}

	// check for other scheduler-related args
	for(++a; a < argc; ++a)
	{
		// gpu
		if(!strcmp(argv[a], "--yield"))
		{
			m_gpuYield = true;
		}
		else if(!strcmp(argv[a], "--anaGPU"))
		{
			m_fAnalysisOnGPU = true;
		}
		else if(!strcmp(argv[a], "--byRank"))
		{
			if(!getArg(m_setGpuByMpiRank, ++a, argc, argv))
			{
				m_setGpuByMpiRank = 0;
				--a;
			}
		}
		else if(!strcmp(argv[a], "-l") || !strcmp(argv[a], "--localLayout"))
		{
			if(!getArg(m_localLayoutPreset, ++a, argc, argv))
			{
				std::cerr << "[EE][SchedulerConfig] expected local layout preset id for --localLayout\n";
				return -1;
			}
			if(m_localLayoutPreset == "dis");
			else if(m_localLayoutPreset == "1_1T4_4");
			else if(m_localLayoutPreset == "2_2");
			else if(m_localLayoutPreset == "2_2T4_3");
			else if(m_localLayoutPreset == "2_2T4_4");
			else if(m_localLayoutPreset == "4_1T4_3");
			else if(m_localLayoutPreset == "3_1T4_3");
			else if(m_localLayoutPreset == "1_3T4_3");
			else if(m_localLayoutPreset == "3_2T4_3");
			else if(m_localLayoutPreset == "default");
			else if(m_localLayoutPreset == "dis_div2");
			else
			{
				std::cerr << "[EE][SchedulerConfig] no such local layout preset: " << m_localLayoutPreset << '\n';
				return -1;
			}
		}

		// cpu
		else if(!strcmp(argv[a], "--1d"))
		{
			m_dimensions = 1;
		}
		else if(!strcmp(argv[a], "--mcsThreads"))
		{
			if(!getArg(m_mcsThreads, ++a, argc, argv))
			{
				std::cerr << "Expected number of threads for --mcsThreads.\n";
				return -1;
			}
			if(m_mcsThreads < 1)
			{
				std::cerr << "Invalid number of threads given, irgnoring.\n";
				m_mcsThreads = 1;
			}
		}
		else if(!strcmp(argv[a], "--ddBlocksize"))
		{
			if(!getArg(m_ddBlocksize[0], ++a, argc, argv))
			{
				std::cerr << "Expected ddBlocksize for --ddBlocksize.\n";
				return -1;
			}
			if(!getArg(m_ddBlocksize[1], ++a, argc, argv))
			{
				--a;
				m_ddBlocksize[1] = m_ddBlocksize[0];
			}
		}

		else
		{
			--a;
			break;
		}
	}
	if(a >= argc)
		return --a;

	return a;
}

Kpz::SchedulerService* Kpz::SchedulerConfig::mkScheduler(int x, int y) const
{
	int preferredSchedulerType = m_preferredSchedulerType;
	if(preferredSchedulerType == SCHEDULER_TYPE_ANY)
	{
#ifdef MODE_CUDA
		preferredSchedulerType = SCHEDULER_TYPE_GPU;
#else
		preferredSchedulerType = SCHEDULER_TYPE_CPU;
#endif
	}

	if(preferredSchedulerType == SCHEDULER_TYPE_CPU)
		return mkSchedulerCPU(x,y);
	else
		return mkSchedulerGPU(x,y);
}

Kpz::SchedulerService* Kpz::SchedulerConfig::mkScheduler(Kpz::SchedulerService&& other) const
{
	int preferredSchedulerType = m_preferredSchedulerType;
	if(preferredSchedulerType == SCHEDULER_TYPE_ANY)
	{
#ifdef MODE_CUDA
		preferredSchedulerType = SCHEDULER_TYPE_GPU;
#else
		preferredSchedulerType = SCHEDULER_TYPE_CPU;
#endif
	}

	if(preferredSchedulerType == SCHEDULER_TYPE_CPU)
	{
		return mkSchedulerCPU(std::move(other));
	}
	else
		return mkSchedulerGPU(std::move(other));
}

#include "kpzSchedulersCPU.h"
Kpz::SchedulerService* Kpz::SchedulerConfig::mkSchedulerCPU(int x, int y) const
{
	SchedulerService* scheduler = 0;
	if( m_dimensions == 2)
	{
		if(m_schedulerId == "db")
		{
			auto tmpScheduler = new Kpz::SchedulerCPUDB<Kpz::RngIface::DSFMT>(x, y, m_mcsThreads);
			scheduler = tmpScheduler;
			if(m_db_lORest >= 0)
			{
				tmpScheduler->setlORest(m_db_lORest);
			}
			if(m_ddBlocksize[0] > 0 && m_ddBlocksize[1] > 0)
			{
				tmpScheduler->setLayout(m_ddBlocksize);
			}
			tmpScheduler->printLayout(std::cout);
		}
		else if(m_schedulerId == "dtStatic")
		{
			auto tmpScheduler = new Kpz::SchedulerCPUDT<Kpz::RngIface::DSFMT>(x, y, m_mcsThreads);
			scheduler = tmpScheduler;
			if(m_ddBlocksize[0] > 0 && m_ddBlocksize[1] > 0)
			{
				tmpScheduler->setLayout(m_ddBlocksize);
			}
			tmpScheduler->printLayout(std::cout);
			if(!tmpScheduler->checkLayout())
			{
				std::cerr << "[EE][mkSchedulerCPU] SchedulerCPUDT: invalid DD layout.\n";
				delete tmpScheduler;
				scheduler = 0;
			}
		}
		else if(m_schedulerId == "dt")
		{
			auto tmpScheduler = new Kpz::SchedulerCPUDTr<Kpz::RngIface::TinyMT32>(x, y, m_mcsThreads);
			scheduler = tmpScheduler;
			if(m_ddBlocksize[0] > 0 && m_ddBlocksize[1] > 0)
			{
				tmpScheduler->setLayout(m_ddBlocksize);
			}
			tmpScheduler->printLayout(std::cout);
			if(!tmpScheduler->checkLayout())
			{
				std::cerr << "[EE][mkSchedulerCPU] SchedulerCPUDTr: invalid DD layout.\n";
				delete tmpScheduler;
				scheduler = 0;
			}
		}
		else if(m_schedulerId == "sca")
		{
			scheduler = new Kpz::SchedulerCPUCB<Kpz::RngIface::TinyMT32>(x, y, m_mcsThreads);
		}
		else if(m_schedulerId == "scaBit")
		{
			scheduler = new Kpz::SchedulerCPUCBBit<Kpz::RngIface::TinyMT32>(x, y, m_mcsThreads);
		}
		else
			scheduler = new Kpz::SchedulerCPU(x, y);
	}
	else if (m_dimensions == 1)
		scheduler = new Kpz::SchedulerCPU1D(x, y);
	else
	{
		std::cerr << "BUG: dimension invalid.\n";
	}
	setFlags(scheduler);
	return scheduler;
}

Kpz::SchedulerService* Kpz::SchedulerConfig::mkSchedulerCPU(SchedulerService&& other) const
{
	std::cerr << "[BUG][SchedulerConfig] move verison of mkSchedulerCPU not implemented.\n";
	return 0;
}

#ifdef MODE_CUDA
#include "kpzSchedulersGPU.h"
#endif

Kpz::SchedulerService* Kpz::SchedulerConfig::mkSchedulerGPU(int x, int y) const
{
#ifdef MODE_CUDA
	const int device = initCUDA();
	SchedulerBase* scheduler = 0;
	if(m_schedulerId == "sca")
	{
		scheduler = new Kpz::SchedulerSCA< KPZ_SWITCH_RNG >(x, y, device);
	}
	else if(m_schedulerId == "scaBit")
	{
		scheduler = new Kpz::SchedulerSCABit<KPZ_SWITCH_RNG>(x, y, device);
	}
	else if(m_schedulerId == "db4")
	{
		scheduler = new Kpz::Scheduler< KPZ_SWITCH_RNG, Scheduler_LocalLayoutDefault >(x, y, device);
	}
	else
	{
		if(m_localLayoutPreset == "dis")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutDis >(x, y, device);
		else if(m_localLayoutPreset == "2_2")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutTC2_2 >(x, y, device);
		else if(m_localLayoutPreset == "2_2T4_3")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutTC2_2T4_3 >(x, y, device);
		else if(m_localLayoutPreset == "2_2T4_4")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutTC2_2T4_4 >(x, y, device);
		else if(m_localLayoutPreset == "4_1T4_3")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutTC4_1T4_3 >(x, y, device);
		else if(m_localLayoutPreset == "3_1T4_3")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutTC3_1T4_3 >(x, y, device);
		else if(m_localLayoutPreset == "1_3T4_3")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutTC1_3T4_3 >(x, y, device);
		else if(m_localLayoutPreset == "3_2T4_3")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutTC3_2T4_3 >(x, y, device);
		else if(m_localLayoutPreset == "dis_div2")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutDis_div2 >(x, y, device);
		else if(m_localLayoutPreset == "1_1T4_4")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutTC1_1T4_4 >(x, y, device);
		else
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutDefault >(x, y, device);
	}
	setFlags(scheduler);
	return scheduler;
#else
	std::cerr << "Cannot create GPU scheduler in this version!\n";
	exit(1);
#endif
}

Kpz::SchedulerService* Kpz::SchedulerConfig::mkSchedulerGPU(SchedulerService&& other) const
{
#ifdef MODE_CUDA
	const int device = initCUDA();
	SchedulerBase* scheduler = 0;
	if(m_schedulerId == "sca")
	{
		scheduler = new Kpz::SchedulerSCA< KPZ_SWITCH_RNG >(std::move(other), device);
	}
	else if(m_schedulerId == "db4")
	{
		scheduler = new Kpz::Scheduler< KPZ_SWITCH_RNG, Scheduler_LocalLayoutDefault >(std::move(other), device);
	}
	else
	{
		if(m_localLayoutPreset == "dis")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutDis >(std::move(other), device);
		else if(m_localLayoutPreset == "2_2")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutTC2_2 >(std::move(other), device);
		else if(m_localLayoutPreset == "2_2T4_3")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutTC2_2T4_3 >(std::move(other), device);
		else if(m_localLayoutPreset == "2_2T4_4")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutTC2_2T4_4 >(std::move(other), device);
		else if(m_localLayoutPreset == "4_1T4_3")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutTC4_1T4_3 >(std::move(other), device);
		else if(m_localLayoutPreset == "3_1T4_3")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutTC3_1T4_3 >(std::move(other), device);
		else if(m_localLayoutPreset == "1_3T4_3")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutTC1_3T4_3 >(std::move(other), device);
		else if(m_localLayoutPreset == "3_2T4_3")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutTC3_2T4_3 >(std::move(other), device);
		else if(m_localLayoutPreset == "dis_div2")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutDis_div2 >(std::move(other), device);
		else if(m_localLayoutPreset == "1_1T4_4")
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutTC1_1T4_4 >(std::move(other), device);
		else
			scheduler = new Kpz::SchedulerDT< KPZ_SWITCH_RNG, SchedulerDT_LocalLayoutDefault >(std::move(other), device);
	}
	setFlags(scheduler);
	return scheduler;
#else
	std::cerr << "Cannot create GPU scheduler in this version!\n";
	exit(1);
#endif
}

const std::string& Kpz::SchedulerConfig::schedulerId() const
{
#ifdef MODE_CUDA
	return m_schedulerId;
#else
	return m_schedulerId;
#endif
}

Kpz::FrontendArgs::FrontendArgs(int nSchedulerConfigs)
	: m_h5in(0), m_h5out(0)
	, m_fSave(false)
	, m_fForcePrevMeasureOverride(false)
	, m_fCountTagsAsMeasurement(false)
	, m_sSamplingSuspended(false), m_sFinished(false)
	, m_disorderFill(NAN)
	, m_schedulerConfig(nSchedulerConfigs)
{
	for(int a = 0; a < 4; ++a)
		m_rateQ[a] = m_rateP[a] = NAN;
}

int Kpz::FrontendArgs::parseSysArgs(int argc, char* argv[])
{
	int nextArg = 1;
	if(!getArg(m_mcs, nextArg, argc, argv))
	{
		std::cerr << "Number of mcs was not supplied.\n";
		return -1;
	}
	if(!getArg(m_x, ++nextArg, argc, argv))
	{
		m_x = m_y = 10; // set default system size
		return nextArg;
	}
	if(!getArg(m_y, ++nextArg, argc, argv))
	{
		m_y = m_x;
		return nextArg;
	}
	return nextArg+1;
}

int Kpz::FrontendArgs::parseArg(int argc, char* argv[], int a)
{
	if((!strcmp(argv[a], "--scheduler")) || (!strcmp(argv[a], "--gpuscheduler")) || (!strcmp(argv[a], "--cpuscheduler")))
	{
		int preferredSchedulerType = SchedulerConfig::SCHEDULER_TYPE_ANY;
		switch(argv[a][2])
		{
			case 'g':
				preferredSchedulerType = SchedulerConfig::SCHEDULER_TYPE_GPU;
				break;
			case 'c':
				preferredSchedulerType = SchedulerConfig::SCHEDULER_TYPE_CPU;
				break;
		}
		int tmp;
		if(!getArg(tmp, ++a, argc, argv))
		{
			tmp = 0;
			--a;
		}
		return m_schedulerConfig.at(tmp).parseArgs(argc,argv,++a, preferredSchedulerType);
	}
	else if(!strcmp(argv[a], "--nschedulers"))
	{
		int tmp;
		if(!getArg(tmp, ++a, argc, argv))
		{
			std::cerr << "[EE][FrontendArgs] expected number of schedulers for --nschedulers.\n";
			return -1;
		}
		if(tmp < m_schedulerConfig.size())
			std::cerr << "[WW][FrontendArgs] Potentiall discarding scheduler configs.\n";
		m_schedulerConfig.resize(tmp);
	}
	else if(!strcmp(argv[a], "--disorder2"))
	{
		if(!getArg(m_disorderFill, ++a, argc, argv))
		{
			m_disorderFill = 0.5;
			--a;
		}
	}
	else if(!memcmp(argv[a], "--pq", 4))
	{
		int n = 0;
		if(strlen(argv[a]) > 4)
		{
			char tmp[2];
			tmp[1] = 0;
			tmp[0] = argv[a][4];
			n = atoi(tmp);
		}
		if(!getArg(m_rateP[n], ++a, argc, argv))
		{
			std::cerr << "--pq: Expected a double value for p.\n";
			return -1;
		}
		if(!getArg(m_rateQ[n], ++a, argc, argv))
		{
			m_rateQ[n] = 0.;
			--a;
		}
	}
	else if(!strcmp(argv[a], "--anaThreads"))
	{
		if(!getArg(m_anaThreads, ++a, argc, argv))
		{
			std::cerr << "Expected number of threads for --anaThreads.\n";
			return -1;
		}
		if(m_anaThreads < 1)
		{
			std::cerr << "Invalid number of threads given, irgnoring. (--anaThreads)\n";
			m_anaThreads = 1;
		}
	}
	else if(!strcmp(argv[a], "--mcsFct"))
	{
		const char* tmp;
		if(!getArg(tmp, ++a, argc, argv))
		{
			std::cerr << "[EE][FrontendArgs] Expected an idetifier for --mcsFct.\n";
			return -1;
		}
		if(!strcmp(tmp, "nop"))
		{
			m_mcsFct = MCSFCT_NOP;
		}
		else
		{
			std::cerr << "[EE][FrontendArgs] Unknown id for --mcsFct: " << tmp << '\n';
			return -1;
		}
	}
	else if(!strcmp(argv[a], "--prevMeasureOverride"))
	{
		unsigned int ovTime;
		if(!getArg(ovTime, ++a, argc, argv))
		{
			std::cerr << "Expected an integer value for --prevMeasureOverride.\n";
			return -1;
		}
		m_prevMeasureOverrideTime.push(ovTime);
		const char* tmp;
		if(getArg(tmp, ++a, argc, argv) && (!strcmp(tmp, "-f")))
			m_fForcePrevMeasureOverride = true;
		else
			--a;
	}
	else if(!strcmp(argv[a], "--noSave"))
	{
		m_fSave = false;
	}
	else if(!strcmp(argv[a], "--save"))
	{
		m_fSave = true;
	}
	else if(!strcmp(argv[a], "--countTagsAsMeasurement"))
	{
		m_fCountTagsAsMeasurement = true;
	}
	else if(!strcmp(argv[a], "--noCountTagsAsMeasurement"))
	{
		m_fCountTagsAsMeasurement = false; // default
	}
	else if(!strcmp(argv[a], "--saveH5"))
	{
		if(!getArg(m_h5out, ++a, argc, argv))
		{
			std::cerr << "Expected a HDF5 filname-prefix for --saveH5.\n";
			return -1;
		}
	}
	else if(!strcmp(argv[a], "--loadH5"))
	{
		if(!getArg(m_h5in, ++a, argc, argv))
		{
			std::cerr << "Expected a HDF5 filname-prefix for --loadH5.\n";
			return -1;
		}
	}
	else if(!strcmp(argv[a], "--retryH5"))
	{
		if(!getArg(m_h5RetryCount, ++a, argc, argv))
		{
			std::cerr << "Expected a postive integer for --retryH5.\n";
			return -1;
		}
	}

	else if(!strcmp(argv[a], "--maxH"))
	{
		float maxH;
		if((!getArg(maxH, ++a, argc, argv)) || maxH <= 0.)
		{
			std::cerr << "Expected a positive floating point value for --maxH.\n";
			return -1;
		}
		Kmc::TimerSingleton::create(maxH);
	}
	else if(!strcmp(argv[a], "--version"))
	{
		version();
	}

	// sampling
	else if(!strcmp(argv[a], "--sampleRate"))
	{
		double tmp;
		if(!getArg(tmp, ++a, argc, argv))
		{
			std::cerr << "Expected a floating point value for --sampleRate [rate] [startTime=0].\n";
			return 1;
		}
		unsigned int start;
		if(!getArg(start, ++a, argc, argv))
		{
			--a;
			start = 0;
		}
		auto s = SamplingParam::makeExp(start, tmp);
		if(start == 0)
			m_currentSampling = s;
		else
			m_events.push_back(s);
	}
	else if(!strcmp(argv[a], "--linSample"))
	{
		unsigned int tmp;
		if(!getArg(tmp, ++a, argc, argv))
		{
			std::cerr << "Expected an unsigned int value as first arg for --linSample.\n";
			return 1;
		}
		unsigned int start;
		if(!getArg(start, ++a, argc, argv))
		{
			--a;
			start = 0;
		}
		auto s = SamplingParam::makeLinear(start, tmp);
		if(start == 0)
			m_currentSampling = s;
		else
			m_events.push_back(s);
	}
	else if(!strcmp(argv[a], "--randomSample"))
	{
		double tmp;
		if(!getArg(tmp, ++a, argc, argv))
		{
			std::cerr << "Expected a floating point value for --randomSample.\n";
			return 1;
		}
		unsigned int start;
		if(!getArg(start, ++a, argc, argv))
		{
			--a;
			start = 0;
		}
		auto s = SamplingParam::makeRandom(start, tmp/2.); // uniform distribution [mean/2:mean+mean/2]
		if(start == 0)
			m_currentSampling = s;
		else
			m_events.push_back(s);
	}
	else if(!strcmp(argv[a], "--corrStart"))
	{
		unsigned int t, s;
		if(!getArg(t, ++a, argc, argv))
		{
			std::cerr << "[EE][FrontendArgs] Expected a unsigned int for --corrStart.\n";
			return 1;
		}
		if(!getArg(s, ++a, argc, argv))
		{
			--a;
			s = std::numeric_limits<unsigned int>::max();
		}
		m_events.push_back(std::shared_ptr<SamplingEvent>(new CorrelateTag(t,s, fCountTagsAsMeasurement())));
	}

	else
	{
		std::cerr << "[EE][FrontendArgs] Unknown argument: " << argv[a] << '\n';
		return -1;
	}
	return a;
}

void Kpz::FrontendArgs::init()
{
	std::cout << "[MM][FrontendArgs] running " << m_anaThreads << " threads for analysis (per process)\n";
	std::cout << std::setprecision(12);
	std::cout << "[MM][FrontendArgs] output precision= " << 12 << '\n'; // TODO: option
	Kpz::Correlator::setNThreads(m_anaThreads);
	m_events.push_back(std::shared_ptr<SamplingEvent>(new EndEvent(m_mcs, false)));
	m_events.sort([](const std::shared_ptr<SamplingEvent>& a, const std::shared_ptr<SamplingEvent>& b){return a->startTime() < b->startTime();});
	if(!m_currentSampling)
		m_currentSampling = SamplingParam::makeExp(0, .1);
}

bool Kpz::FrontendArgs::overridePrevMeasurementTime(SchedulerService* scheduler)
{
	fastForwardEvents(scheduler->mcs());
	for(; prevMeasureOverride(); popPrevMeasurementOverride())
	{
		unsigned int nextEval;
		while(prevMeasureOverride() && prevMeasureOverrideTime() > scheduler->mcs())
			popPrevMeasurementOverride();
		if(!prevMeasureOverride())
			return false;
		if(m_fForcePrevMeasureOverride)
		{
			nextEval = sampling().nextEvalFastForward(prevMeasureOverrideTime(), scheduler);
		}
		else
		{
			nextEval = sampling().nextEval(prevMeasureOverrideTime(), scheduler);
			if(nextEval < scheduler->mcs())
			{
				std::cerr << "[WW] prev measure-time override leads to next measurement being in the past: ignoring.\n";
				continue;
			}
		}
		for(auto a = m_events.begin(); a != m_events.end(); ++a)
		{
			if((*a)->startTime() > nextEval)
			{
				// insert placeholder event, so that a measurement will occour at this time
				m_events.insert(a, std::shared_ptr<SamplingEvent>(new SamplingBarrier(nextEval, true)));
				break;
			}
			else if((*a)->fCountAsMeasurement())
			{
				std::cerr << "[WW] prev measure-time override overridden by intermediate sampling event.\n";
				return false;
			}
		}
		setSamplingSuspended();
		return true;
	}
	return false;
}

void Kpz::FrontendArgs::fastForwardEvents(unsigned int nowMcs)
{
	while((!m_events.empty()) && m_events.front()->startTime() < nowMcs)
	{
		auto s = std::dynamic_pointer_cast<SamplingParam>(m_events.front());
		if(s)
			m_currentSampling = s;
		m_events.pop_front();
	}
}

bool Kpz::FrontendArgs::handleNextEvent(SchedulerService* scheduler)
{
	if(m_events.empty())
		return true;
	auto e = m_events.front();
	bool countTagAsMeasurement = (e->startTime() > scheduler->mcs());
	while(e->startTime() <= scheduler->mcs())
	{
		m_events.pop_front();
		if(e->startTime() == scheduler->mcs())
		{
			e->handle(*this, scheduler);
			countTagAsMeasurement |= e->fCountAsMeasurement();

		}
		else
			std::cerr << "[WW][FrontendArgs][handleNextEvent] Event skipped: "
				<< e->startTime() << '<' << scheduler->mcs() << '\n';
		if(m_events.empty())
			break;
		else
			e = m_events.front();
	}
	return countTagAsMeasurement;
}

void Kpz::FrontendArgs::insertEvent(std::shared_ptr<SamplingEvent> event)
{
	auto a = std::find_if(m_events.begin(), m_events.end()
			, [&event](const std::shared_ptr<SamplingEvent>& a){return a->startTime() > event->startTime();});
	m_events.insert(a, event);
}

unsigned int Kpz::FrontendArgs::nextMcs(unsigned int mcs, SchedulerService* scheduler) const
{
	unsigned int ret = sampling().nextEval(mcs, scheduler);
	if(!m_events.empty())
	{
		ret = std::min(m_events.front()->startTime(), ret);
	}
	return ret - scheduler->mcs();
}

bool Kpz::FrontendArgs::setPQ(SchedulerService* scheduler) const
{
	for(int a = 0; a < 4 && !std::isnan(m_rateP[a]); ++a)
		if(!scheduler->setPQ(m_rateP[a], m_rateQ[a], a))
			return false;
	return true;
}
