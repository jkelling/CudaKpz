/***************************************************************************
*   Copyright 2011 - 2012 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef KPZ_SCHEDULER_CPU_CELL_H
#define KPZ_SCHEDULER_CPU_CELL_H

#include "schedulerServiceStash.h"

#include "LocalLayout.h"

namespace Kpz
{

	template<class LocalLayout>
	class SchedulerCPUCell : public SchedulerServiceStash
	{
		protected:

		void updateThreadCell(int blocks, int xw, int yw);
		void updateThreadCell(int* system);

		void init();

		public:
		typedef LocalLayoutProxy<LocalLayout, LocalLayoutDB4> MyLocalLayout;
		static const int ACT_THREAD_CELL_DIM_X = MyLocalLayout::THREAD_CELL_DIM_X - DEAD_BORDER_DIM;
		static const int ACT_THREAD_CELL_DIM_Y = MyLocalLayout::THREAD_CELL_DIM_X - DEAD_BORDER_DIM;

			SchedulerCPUCell(int lx, int ly) : SchedulerServiceStash(lx, ly) {
				init();
			}
			SchedulerCPUCell(const SystemSize& size) : SchedulerServiceStash(size) {
				init();
			}

		void mcs(int n);

	};

}

#endif
