/***************************************************************************
*   Copyright 2011 - 2014 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef KPZ_SCHEDULER_H
#define KPZ_SCHEDULER_H

#include "schedulerBase.h"

#include "kpzConst.h"
#include "LocalLayout.h"

namespace Kpz
{
#ifdef KPZ_FERMI
		typedef LocalLayout<5, 5, 1, 2, 1> Scheduler_LocalLayoutDefault;
#else
		typedef LocalLayout<5, 4, 1, 1, 1> Scheduler_LocalLayoutDefault;
#endif

	template<class GpuRng, class LocalLayout>
	class Scheduler : public SchedulerBase
	{
		protected:
		GpuRng m_gpuRng;
		int m_randomNumberBlockMultiplier;
		int calcRandomNumberBlockMultiplier() const {
			int t = m_size.blocksGPU()/m_blocks;
			if(m_size.blocksGPU()%m_blocks > 0)
				++t;
			return t;
		}

		typedef void (Scheduler::*kernelCallerFct) (int xw, int yw);
		void caller_mcsSyncRng(int xw, int yw);
		void caller_mcsPQSyncRng(int xw, int yw);
		// void caller_mcsDisorder2SyncRng(int xw, int yw);

		void doMcs(int n, kernelCallerFct mcs);

		void init();
		virtual bool changedSizeEvent();

		public:
		typedef LocalLayoutProxy<LocalLayout, LocalLayoutDB4> MyLocalLayout;

			Scheduler(int lx, int ly, int device = -1);
			Scheduler(const SystemSize& size, int device = -1);
			Scheduler(SchedulerService &&other, int device = -1);

		void randomize() {m_gpuRng.randomize();}

		void mcs(int n) {
			doMcs(n, &Scheduler::caller_mcsSyncRng);
		}
		void mcsSyncRng(int n) {
			doMcs(n, &Scheduler::caller_mcsSyncRng);
		}

		void mcsPQSyncRng(int n) {
			doMcs(n, &Scheduler::caller_mcsPQSyncRng);
		}
		void mcsDisorder2SyncRng(int n) {
			// doMcs(n, caller_mcsDisorder2SyncRng);
			std::cerr << "[BUG] mcsDisorder2SyncRng() not implemented.\n";
		}

		virtual double mcsScale() const;

		virtual bool collectH5(splash::DataCollector* data, const std::string& prefix = "");
		virtual bool retrieveH5(splash::DataCollector* data, const std::string& prefix = "");
	};
}
#endif
