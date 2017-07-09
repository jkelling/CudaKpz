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

#include "GPURandom123_Threefry.h"

#define KMC_L_LCG_N 4
#include "../kmcRandom.h"

#include <cudaError.h>

__shared__ unsigned int* Random123_Threefry_Device_d_Random;
__constant__ Random123_Threefry::Threefry::key_type Random123_Threefry_Device_Key;

/*! \warning, an __syncthreads(); should occour sometime between constructoring
 * and destructing for the storing of state to work: It is required, that thread(0) rand before any other thread hits ~Device().
 */
class Random123_Threefry::Device
{
	Threefry::ctr_type m_ctr;
	Threefry m_rng;

	public:

		__device__ Device(void* d) {
		m_ctr[0] = threadIdx.x + ((blockIdx.x)*blockDim.x);
		if(threadIdx.x == 0)
			Random123_Threefry_Device_d_Random = (unsigned int*)d;
		m_ctr.v[1] = ((int*)d)[m_ctr[0]];
	}

		__device__ ~Device() {
		Random123_Threefry_Device_d_Random[m_ctr[0]] = m_ctr[1];
	}

	__device__ int random() {
		++m_ctr[1];
		return m_rng(m_ctr, Random123_Threefry_Device_Key)[0];
	}

	__device__ float randomFloat() {
		return (random()&KMC_M_LCG_RAND_RED)/(float)KMC_LCG_RAND_SUP_RED;
	}
};

inline void Random123_Threefry::initDeviceConst()
{
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(Random123_Threefry_Device_Key, &m_key, sizeof(m_key)));
}
