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

#ifndef KPZ_KPZ_RNG
#define KPZ_KPZ_RNG

#include "GPURandom.h"
#include "GPUTinyMT.h"

#include "kpzConst.h"
#include "cudaDefines.h"

#include <tinyMT/tinymt32_kernel.cuh>

class SLCG64::Device
{
	unsigned long long i;
	unsigned long long* d_Random;

	public:

	//! \todo coalescing
		__device__ Device(void* d) : d_Random((unsigned long long*)d) {
		const int THREAD_RNG_OFFSET =
			(GPU_BLOCK_ID_X *( GPU_BLOCK_DIM_Y ? GPU_BLOCK_DIM_X*GPU_BLOCK_DIM_Y : GPU_BLOCK_DIM_X)) \
			+ GPU_THREAD_ID_X + ((GPU_THREAD_ID_Y)*GPU_BLOCK_DIM_X);
		i = d_Random[THREAD_RNG_OFFSET];
	}

		__device__ ~Device() {
		const int THREAD_RNG_OFFSET =
			(GPU_BLOCK_ID_X *( GPU_BLOCK_DIM_Y ? GPU_BLOCK_DIM_X*GPU_BLOCK_DIM_Y : GPU_BLOCK_DIM_X)) \
			+ GPU_THREAD_ID_X + ((GPU_THREAD_ID_Y)*GPU_BLOCK_DIM_X);
		d_Random[THREAD_RNG_OFFSET] = i;
	}

	__device__ int random() {
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

class TinyMT32::Device
{
	unsigned int* d_Random;
	TINYMT32_STATUS_T m_state;

	public:

		__device__ Device(void* d) : d_Random((unsigned int*)d) {
		int THREAD_RNG_OFFSET =
			(GPU_BLOCK_ID_X *( GPU_BLOCK_DIM_Y ? GPU_BLOCK_DIM_X*GPU_BLOCK_DIM_Y : GPU_BLOCK_DIM_X)) \
			+ GPU_THREAD_ID_X + ((GPU_THREAD_ID_Y)*GPU_BLOCK_DIM_X);
		const int N_GENERATORS = GPU_GRID_DIM_X*( GPU_BLOCK_DIM_Y ? GPU_BLOCK_DIM_X*GPU_BLOCK_DIM_Y : GPU_BLOCK_DIM_X);
#pragma unroll
		for(int a = 0; a < 4; ++a)
		{
			m_state.status[a] = d_Random[THREAD_RNG_OFFSET];
			THREAD_RNG_OFFSET += N_GENERATORS;
		}
		m_state.mat1 = d_Random[THREAD_RNG_OFFSET];
		THREAD_RNG_OFFSET += N_GENERATORS;
		m_state.mat2 = d_Random[THREAD_RNG_OFFSET];
		THREAD_RNG_OFFSET += N_GENERATORS;
		m_state.tmat = d_Random[THREAD_RNG_OFFSET];
	}

		__device__ ~Device() {
		int THREAD_RNG_OFFSET =
			(GPU_BLOCK_ID_X *( GPU_BLOCK_DIM_Y ? GPU_BLOCK_DIM_X*GPU_BLOCK_DIM_Y : GPU_BLOCK_DIM_X)) \
			+ GPU_THREAD_ID_X + ((GPU_THREAD_ID_Y)*GPU_BLOCK_DIM_X);
		const int N_GENERATORS = GPU_GRID_DIM_X*( GPU_BLOCK_DIM_Y ? GPU_BLOCK_DIM_X*GPU_BLOCK_DIM_Y : GPU_BLOCK_DIM_X);
		for(int a = 0; a < 4; ++a)
		{
			d_Random[THREAD_RNG_OFFSET] = m_state.status[a];
			THREAD_RNG_OFFSET += N_GENERATORS;
		}
	}

	__device__ int random() {
		return (int)tinymt32_uint32(&m_state);
	}

	__device__ float randomFloat() {
		return tinymt32_single(&m_state);
	}

	__device__ void seed() {
		tinymt32_init(&m_state, m_state.status[0]);
	}
};

#endif
