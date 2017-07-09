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

#ifndef KPZ_SCHEDULER_CPU_DT_H
#define KPZ_SCHEDULER_CPU_DT_H

#include "schedulerCPUpThread.h"
#include "systemSize.h"

#include <doubleTilingDD.h>

namespace Kpz
{

	template<class Rng>
	class SchedulerCPUDT : public SchedulerCPUpThreadRng<Rng>
	{
		protected:
		using SchedulerCPU::m_deviceMcs;
		using SchedulerCPU::m_dsfmt;
		using SchedulerCPUpThread::run;

		typedef Kmc::DoubleTilingDD<2> DD;
		DD m_dd;
		
		void layoutFromL1();

		static const int L_MIN[2];

		struct Func_MCSBase;

		public:
		
			SchedulerCPUDT(int lx, int ly, int threads = 1);
			SchedulerCPUDT(const SystemSize& size, int threads = 1);

		virtual void mcs(int n);
		// virtual void mcsPQ(int n);
		// virtual void mcsDisorder2(int n);

		// virtual void mcsSyncRng(int n);
		// virtual void mcsPQSyncRng(int n);
		// virtual void mcsDisorder2SyncRng(int n);
		
		void setLayout(const int layout[2]);
		std::ostream& printLayout(std::ostream& out = std::cout);
		bool checkLayout() const {return m_dd.nBlocks() > 0;}

		using SchedulerCPUpThread::numThreads;
		using SchedulerCPUpThread::joinMCS;
	};

}

#endif
