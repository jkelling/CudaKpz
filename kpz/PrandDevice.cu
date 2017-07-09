/***************************************************************************
*   Copyright 2014 - 2014 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef KPZ_PRAND_DEVICE
#define KPZ_PRAND_DEVICE

#include "Prand.h"

#include "kpzConst.h"
#include "cudaDefines.h"

class PrandMT::Device
{
	unsigned int i;
	unsigned int* d_Random;

	public:

		__device__ Device(void* d) : d_Random((unsigned int*)d) {
		i = (GPU_BLOCK_ID_X*(( GPU_BLOCK_DIM_Y ? GPU_BLOCK_DIM_X*GPU_BLOCK_DIM_Y : GPU_BLOCK_DIM_X)<<1)) \
			+ GPU_THREAD_ID_X + ((GPU_THREAD_ID_Y) * GPU_BLOCK_DIM_X);
	}

		__device__ ~Device() {
	}

	__device__ int random() {
		const int ret = d_Random[i];
		i += (GPU_GRID_DIM_X*(( GPU_BLOCK_DIM_Y ? GPU_BLOCK_DIM_X*GPU_BLOCK_DIM_Y : GPU_BLOCK_DIM_X)<<1));
		return ret;
	}

	__device__ float randomFloat() {
		return (random()&KMC_M_LCG_RAND_RED)/(float)KMC_LCG_RAND_SUP_RED;
	}
};

#endif
