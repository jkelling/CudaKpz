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

#ifndef KPZ_SCHEDULER_CPU_H
#define KPZ_SCHEDULER_CPU_H

#include "schedulerCPUBase.h"

namespace Kpz
{

	class SchedulerCPU : public SchedulerCPUBase
	{
		protected:

		void updateThreadCell(int blocks, int xw, int yw);
		void updateThreadCell(int* system);
		
		inline void pqOneUpdate(int x, int y);
		inline void pqOneOneUpdate(int x, int y);
		inline void pqUpdate(int x, int y);

		template<class Pred>
		void mcs(int n, Pred pred);

		public:
		
		using SchedulerCPUBase::SchedulerCPUBase;

		void mcs(int n);
		virtual void mcsPQ(int n);
		virtual void mcsDisorder2(int n);

		virtual void mcsSyncRng(int n);
		virtual void mcsPQSyncRng(int n);
		virtual void mcsDisorder2SyncRng(int n);
	};

}

#include "schedulerCPU.inl"

#endif
