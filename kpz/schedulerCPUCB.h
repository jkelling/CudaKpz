/***************************************************************************
*   Copyright 2014 - 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef KPZ_SCHEDULER_CPU_CB_H
#define KPZ_SCHEDULER_CPU_CB_H

#include "schedulerCPUBaseThreadedRng.h"
#include "systemSize.h"

namespace Kpz
{
	/*! CPU scheduler with lattce checkerboard DD (SCA)
	 *
	 * No stash. Will join analysis before beginning MCS.
	 *
	 */
	template<class Rng>
	class SchedulerCPUCB : public SchedulerCPUBaseThreadedRng<Rng>
	{
		protected:

		using SchedulerService::m_deviceMcs;
		using SchedulerService::m_dsfmt;
		using SchedulerService::m_size;
		using SchedulerService::joinRoughness;
		using SchedulerCPUBaseThreaded::run;
		using SchedulerCPUBaseThreaded::joinWorker;
		using SchedulerCPUBaseThreaded::joinAllWorkers;
		using SchedulerCPUBaseThreaded::nextThread;

		int m_signum;

		template<class UpdateFkt>
		struct Func_MCS;
		typedef typename SchedulerCPUBaseThreadedRng<Rng>::template Func_MCSBase<Kpz::SchedulerCPUCB<Rng> > MyFunc_MCSBase;

		template<class UpdateFkt>
		void doMcs(int n, UpdateFkt updateFkt);

		public:
		
			SchedulerCPUCB(int lx, int ly, int threads = 1);
			SchedulerCPUCB(const SystemSize& size, int threads = 1);

		virtual void mcs(int n);

		virtual void mcsSyncRng(int n);
		virtual void mcsPQSyncRng(int n);
		virtual void mcsDisorder2SyncRng(int n);

		using SchedulerCPUBaseThreaded::nThreads;
	};

}

#endif
