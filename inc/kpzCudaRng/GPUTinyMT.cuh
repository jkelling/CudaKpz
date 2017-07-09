/***************************************************************************
*   Copyright 2014 - 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#pragma once

#include "GPUTinyMT.h"

#include "../tinyMT/tinymt32_kernel.cuh"

__shared__ unsigned int* TinyMT32_Device_d_Random;

/*! \warning, an __syncthreads(); should occour sometime between constructoring
 * and destructing for the storing of state to work: It is required, that thread(0) rand before any other thread hits ~Device().
 */
class TinyMT32::Device
{
	TINYMT32_STATUS_T m_state;

	public:

		__device__ Device(void* d) {
		int THREAD_RNG_OFFSET = threadIdx.x + ((threadIdx.y)*blockDim.x);
		if(THREAD_RNG_OFFSET == 0)
			TinyMT32_Device_d_Random = (unsigned int*)d;
		THREAD_RNG_OFFSET += (blockIdx.x *( blockDim.y ? blockDim.x*blockDim.y : blockDim.x));
		const int N_GENERATORS = gridDim.x*( blockDim.y ? blockDim.x*blockDim.y : blockDim.x);
#pragma unroll
		for(int a = 0; a < 4; ++a)
		{
			m_state.status[a] = ((unsigned int*)d)[THREAD_RNG_OFFSET];
			THREAD_RNG_OFFSET += N_GENERATORS;
		}
		m_state.mat1 = ((unsigned int*)d)[THREAD_RNG_OFFSET];
		THREAD_RNG_OFFSET += N_GENERATORS;
		m_state.mat2 = ((unsigned int*)d)[THREAD_RNG_OFFSET];
		THREAD_RNG_OFFSET += N_GENERATORS;
		m_state.tmat = ((unsigned int*)d)[THREAD_RNG_OFFSET];
	}

		__device__ ~Device() {
		int THREAD_RNG_OFFSET =
			(blockIdx.x *( blockDim.y ? blockDim.x*blockDim.y : blockDim.x)) \
			+ threadIdx.x + ((threadIdx.y)*blockDim.x);
		const int N_GENERATORS = gridDim.x*( blockDim.y ? blockDim.x*blockDim.y : blockDim.x);
		for(int a = 0; a < 4; ++a)
		{
			TinyMT32_Device_d_Random[THREAD_RNG_OFFSET] = m_state.status[a];
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
