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

#ifndef KPZ_SCHEDULER_CPU_DTR_H
#define KPZ_SCHEDULER_CPU_DTR_H

#include "schedulerCPUBaseThreadedRng.h"

#include <doubleTilingDD.h>

namespace Kpz
{
	/*! CPU scheduler with random-origin DT-DD.
	 *
	 * No stash. Will join analysis before beginning MCS.
	 *
	 */
	template<class Rng>
	class SchedulerCPUDTr : public SchedulerCPUBaseThreadedRng<Rng>
	{
		protected:
		using SchedulerService::m_deviceMcs;
		using SchedulerService::m_dsfmt;
		using SchedulerService::joinRoughness;
		using SchedulerCPUBaseThreaded::run;
		using SchedulerCPUBaseThreaded::joinWorker;
		using SchedulerCPUBaseThreaded::joinAllWorkers;
		using SchedulerCPUBaseThreaded::nextThread;

		typedef Kmc::DoubleTilingDD<2> DD;
		DD m_dd;
		
		void layoutFromL1();

		static const int L_MIN[2];

		template<class UpdateFkt>
		struct Func_MCS;
		typedef typename SchedulerCPUBaseThreadedRng<Rng>::template Func_MCSBase<Kpz::SchedulerCPUDTr<Rng> > MyFunc_MCSBase;

		template<class UpdateFkt>
		void doMcs(int n, UpdateFkt updateFkt);

		public:
		
			SchedulerCPUDTr(int lx, int ly, int threads = 1);
			SchedulerCPUDTr(const SystemSize& size, int threads = 1);

		virtual void mcsSyncRng(int n);
		virtual void mcsPQSyncRng(int n);
		virtual void mcsDisorder2SyncRng(int n);
		
		void setLayout(const int layout[2]);
		std::ostream& printLayout(std::ostream& out = std::cout);
		bool checkLayout() const {return m_dd.nBlocks() > 0;}

		using SchedulerCPUBaseThreaded::nThreads;
	};

}

#endif
