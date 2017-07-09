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

#ifndef KMC_SCHEDULER_GPU
#define KMC_SCHEDULER_GPU

#include "kpzConst.h"

class cudaDeviceProp;

namespace Kpz
{

	class SchedulerGPU
	{
		protected:
			
			int m_blocks, m_device;
			cudaDeviceProp *m_prop;

		inline int getRandomNumberGeneratorCount (int lThreads) {
			return (m_blocks) << lThreads;
		}

		public:
			
			SchedulerGPU (int device = -1, int smemWperBlock = 0);
			~SchedulerGPU ();
	};

}

#endif
