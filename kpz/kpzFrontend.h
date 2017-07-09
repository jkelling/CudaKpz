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

#ifndef KPZ_FRONTEND_H
#define KPZ_FRONTEND_H

#include "schedulerService.h"
#include "sampling.h"

#include <iostream>
#include <cmath>
#include <string>
#include <list>
#include <vector>
#include <queue>
#include <memory>
#include <limits>

namespace Kpz
{
	void version();

	class SchedulerConfig
	{
		public:

		static const short SCHEDULER_TYPE_ANY = 0;
		static const short SCHEDULER_TYPE_GPU = 1;
		static const short SCHEDULER_TYPE_CPU = 2;

		private:
		std::string m_schedulerId;
		int m_db_lORest, m_ddBlocksize[2]; // cpu
		int m_dimensions = 2; // cpu
		int m_mcsThreads = 1; // cpu

		bool m_gpuYield : 1, m_fAnalysisOnGPU : 1; // gpu
		int m_setGpuByMpiRank = -1; // gpu
		std::string m_localLayoutPreset; // gpu

		short m_preferredSchedulerType = SCHEDULER_TYPE_ANY;

		void setFlags(SchedulerService* scheduler) const;

		public:

			SchedulerConfig();

		const std::string& schedulerId() const;
		short preferredSchedulerType() const {return m_preferredSchedulerType;}
		int db_lORest() const {return m_db_lORest;}
		const int* ddBlocksize() const {return m_ddBlocksize;}
		int dimensions() const {return m_dimensions;}
		int mcsThreads() const {return m_mcsThreads;}

		/*!	\return index of next unused arg. -1 on error. */
		int parseArgs(int argc, char* argv[], int a, int preferredSchedulerType = SCHEDULER_TYPE_ANY);

		const std::string& localLayoutPreset() const {return m_localLayoutPreset;}
		void setLocalLayoutPreset(const std::string& llp) {m_localLayoutPreset = llp;}

		SchedulerService* mkScheduler(int x, int y) const;
		SchedulerService* mkScheduler(SchedulerService&& other) const;
		//! \todo
		SchedulerService* mkSchedulerCPU(SchedulerService&& other) const;
		SchedulerService* mkSchedulerCPU(int x, int y) const;

		SchedulerService* mkSchedulerGPU(int x, int y) const;
		SchedulerService* mkSchedulerGPU(SchedulerService&& other) const;

		int initCUDA() const;
	};

	class FrontendArgs
	{
		int m_x, m_y;
		unsigned int m_mcs = MCSFCT_ANY;
		int m_mcsFct = -1;

		double m_rateP[4], m_rateQ[4], m_disorderFill;

		bool m_fSave : 1
			, m_fCountTagsAsMeasurement : 1
			, m_fForcePrevMeasureOverride : 1;
		bool m_sSamplingSuspended : 1
			, m_sFinished : 1;
		const char *m_h5in, *m_h5out;
		unsigned int m_h5RetryCount = 10;
		std::queue<unsigned int> m_prevMeasureOverrideTime;

		std::vector<SchedulerConfig> m_schedulerConfig;

		std::list<std::shared_ptr<SamplingEvent> > m_events;
		std::shared_ptr<SamplingParam> m_currentSampling;
		int m_anaThreads = 1; // cpu threads for parallel analysis

		public:

		FrontendArgs(int nSchedulerConfigs = 1);

		int x() const {return m_x;}
		int y() const {return m_y;}
		unsigned int mcs() const {return m_mcs;}
		int mcsFct() const {return m_mcsFct;}

		double rateP() const {return m_rateP[0];}
		double rateQ() const {return m_rateQ[0];}
		double rateP(int a) const {return m_rateP[a];}
		double rateQ(int a) const {return m_rateQ[a];}
		bool setPQ(SchedulerService* scheduler) const;
		double disorderFill() const {return m_disorderFill;}

		bool fSave() const {return m_fSave;}
		bool fCountTagsAsMeasurement() const {return m_fCountTagsAsMeasurement;}
		bool prevMeasureOverride() const {return m_prevMeasureOverrideTime.size()>0;}
		unsigned int prevMeasureOverrideTime() const {return m_prevMeasureOverrideTime.front();}
		const char* h5in() const {return m_h5in;}
		const char* h5out() const {return m_h5out;}
		unsigned int h5RetryCount() const {return m_h5RetryCount;}

		const SchedulerConfig& schedulerConfig(int n = 0) const {return m_schedulerConfig.at(n);}
		SchedulerConfig& schedulerConfig(int n = 0) {return m_schedulerConfig.at(n);}
		int nSchedulerConfigs() const {return m_schedulerConfig.size();}

		const SamplingParam& sampling() const {return *m_currentSampling;}
		void setCurrentSampling(std::shared_ptr<SamplingParam> sampling) {m_currentSampling = sampling;}
		unsigned int nextMcs(unsigned int mcs, SchedulerService* scheduler) const;
		bool isSamplingSuspended() const {return m_sSamplingSuspended;}
		void setSamplingSuspended(bool b = true) {m_sSamplingSuspended = b;}
		bool isFinished() const {return m_sFinished;}
		void setFinished(bool b = true) {m_sFinished = b;}

		void fastForwardEvents(unsigned int nowMcs);
		unsigned int nextEventTime() const {
			auto n = nextEvent();
			return n ? n->startTime(): std::numeric_limits<unsigned int>::max();
		}
		std::shared_ptr<SamplingEvent> nextEvent() const {
			return (m_events.empty()) ? std::shared_ptr<SamplingEvent>() : m_events.front();
		}
		void popEvent() {if(!m_events.empty()) m_events.pop_front();}
		void insertEvent(std::shared_ptr<SamplingEvent> event);
		/*! \return wether to count this measurement as prevMeasureOverrideTime */
		bool handleNextEvent(SchedulerService* scheduler);
		/*! Set up override for prev measurement time. If override took place,
		 * prevMeasureOverrideTime() will return the used time. Otherwise
		 * prevMeasureOverride() will return false.
		 *
		 * \return wether override took place */
		bool overridePrevMeasurementTime(SchedulerService* scheduler);

		void init();
		void pushPrevMeasureOverrideTime(unsigned int t) {m_prevMeasureOverrideTime.push(t);}
		void popPrevMeasurementOverride() {m_prevMeasureOverrideTime.pop();}
		void clearPrevMeasurementOverride() {m_prevMeasureOverrideTime = std::queue<unsigned int>();}

		/*!	\return index of last used arg. -1 on error. */
		int parseSysArgs(int argc, char* argv[]);
		/*!	\return index of last used arg. -1 on error. */
		int parseArg(int argc, char* argv[], int a);

		static const int MCSFCT_ANY = -1;
		static const int MCSFCT_NOP = 1;
	};

	void registerTimeoutSignals();
	void sigINTtimeout (int sig);
}
#endif
