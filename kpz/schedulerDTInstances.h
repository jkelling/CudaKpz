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

#ifndef KPZ_SCHEDULER_DT_INSTANCES_H
#define KPZ_SCHEDULER_DT_INSTANCES_H

#ifndef KPZ_SWITCH_PRAND
#include "GPURandom.h"
#include "GPUTinyMT.h"

namespace Kpz
{
	template class SchedulerDT<SLCG64, SchedulerDT_LocalLayoutDefault>;
	template class SchedulerDT<TinyMT32, SchedulerDT_LocalLayoutDefault>;
	template class SchedulerDT<SLCG64, SchedulerDT_LocalLayoutDis>;
	template class SchedulerDT<TinyMT32, SchedulerDT_LocalLayoutDis>;
	template class SchedulerDT<TinyMT32, SchedulerDT_LocalLayoutDis_div2>;
	template class SchedulerDT<TinyMT32, SchedulerDT_LocalLayoutTC1_1T4_4>;
	template class SchedulerDT<TinyMT32, SchedulerDT_LocalLayoutTC2_2>;
	template class SchedulerDT<TinyMT32, SchedulerDT_LocalLayoutTC2_2T4_3>;
	template class SchedulerDT<TinyMT32, SchedulerDT_LocalLayoutTC2_2T4_4>;
	template class SchedulerDT<TinyMT32, SchedulerDT_LocalLayoutTC4_1T4_3>;
	template class SchedulerDT<TinyMT32, SchedulerDT_LocalLayoutTC3_1T4_3>;
	template class SchedulerDT<TinyMT32, SchedulerDT_LocalLayoutTC1_3T4_3>;
	template class SchedulerDT<TinyMT32, SchedulerDT_LocalLayoutTC3_2T4_3>;
}
#else
#include "Prand.h"

namespace Kpz
{
	template class SchedulerDT<PrandMT, SchedulerDT_LocalLayoutDefault>;
	template class SchedulerDT<PrandMT, SchedulerDT_LocalLayoutDis>;
}
#endif
#endif
