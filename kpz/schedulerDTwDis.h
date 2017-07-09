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

#pragma once

#include "schedulerDT.h"

namespace Kpz
{
	template<class GpuRng, class LocalLayout>
	class SchedulerDTwDis : public SchedulerDT<GpuRng, LocalLayout>
	{
		typedef SchedulerDT<GpuRng, LocalLayout> Base;

		protected:
		using typename Base::MyLocalLayout;
		using Base::m_gpuRng;
		using Base::m_randomNumberBlockMultiplier;
		using Base::m_disorderP;
		using Base::m_disorderQ;
		using Base::m_blocks;
		using Base::d_System;


		public:

			SchedulerDTwDis(int lx, int ly, int device = -1);
			SchedulerDTwDis(const SystemSize& size, int device = -1);
			SchedulerDTwDis(SchedulerService &&other, int device = -1);

	};
}
