/***************************************************************************
*   Copyright 2011 - 2014 Jeffrey Kelling <j.kelling@hzdr.de>
*                  Helmholtz-Zentrum Dresden-Rossendorf
*                  Institute of Ion Beam Physics and Materials Research
*
*   This program is free software; you can redistribute it and/or
*   modify it under the terms of the GNU General Public
*   License as published by the Free Software Foundation; either
*   version 2 of the License, or (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program; if not, write to the Free Software
*   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
***************************************************************************/

#ifndef KPZ_KPZ_RNG
#define KPZ_KPZ_RNG

#include "GPUSLCG64.h"

#include "../kmcRandom.h"

#include "../tinyMT/tinymt32_kernel.cuh"

class SLCG64::Device
{
	unsigned long long i;
	unsigned long long* d_Random;

	public:

	//! \todo coalescing
		__device__ Device(void* d) : d_Random((unsigned long long*)d) {
		const int THREAD_RNG_OFFSET =
			(blockIdx.x *( blockDim.y ? blockDim.x*blockDim.y : blockDim.x)) \
			+ threadIdx.x + ((threadIdx.y)*blockDim.x);
		i = d_Random[THREAD_RNG_OFFSET];
	}

		__device__ ~Device() {
		const int THREAD_RNG_OFFSET =
			(blockIdx.x *( blockDim.y ? blockDim.x*blockDim.y : blockDim.x)) \
			+ threadIdx.x + ((threadIdx.y)*blockDim.x);
		d_Random[THREAD_RNG_OFFSET] = i;
	}

	static __device__ int random() {
#if 0
		int y = i/KMC_LCG_Q;
		i = KMC_LCG_A*(i - KMC_LCG_Q*y) - KMC_LCG_R*y;
		if(i<0) i+=KMC_LCG_M;
		return i-1;
#else
		i = KMC_SLCG_A_SKIP*i + KMC_SLCG_C_SKIP;
		return (unsigned int)(i^(i>>32))&(-1l); // improve bit randomness avoiding float
#endif
	}

	__device__ float randomFloat() {
		return (random()&KMC_M_LCG_RAND_RED)/(float)KMC_LCG_RAND_SUP_RED;
	}
};
#endif
