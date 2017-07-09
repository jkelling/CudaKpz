/***************************************************************************
*   Copyright 2014 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef KPZ_SCHEDULER_CPU_DB_H
#define KPZ_SCHEDULER_CPU_DB_H

#include "schedulerCPUpThread.h"
#include "systemSize.h"

#include <dynamicOriginDD.h>

namespace Kpz
{

	template<class Rng>
	class SchedulerCPUDB : public SchedulerCPUpThreadRng<Rng>
	{
		protected:
		using SchedulerCPU::m_deviceMcs;
		using SchedulerCPU::m_dsfmt;
		using SchedulerCPUpThread::run;

		typedef Kmc::DynamicOriginDD<2> DD;
		DD m_dd;
		
		void layoutFromL1();

		static const int L_MIN[2];
		int m_lORest[2];

		struct Func_MCSBase;

		public:
		
			SchedulerCPUDB(int lx, int ly, int threads = 1);
			SchedulerCPUDB(const SystemSize& size, int threads = 1);

		virtual void mcs(int n);
		// virtual void mcsPQ(int n);
		// virtual void mcsDisorder2(int n);

		// virtual void mcsSyncRng(int n);
		// virtual void mcsPQSyncRng(int n);
		// virtual void mcsDisorder2SyncRng(int n);
		
		void setlORest(int lr) {m_lORest[0] = m_lORest[1] = lr;}
		void setLayout(const int layout[2]);
		std::ostream& printLayout(std::ostream& out = std::cout);

		using SchedulerCPUpThread::numThreads;
		using SchedulerCPUpThread::joinMCS;
	};

}

#endif
