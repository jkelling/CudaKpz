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

#include "schedulerSCA.h"

#include "kpzConst.h"

#include <kmcMath.h>

#include <iostream>
#include <algorithm>
#include <cassert>

#include <cudaError.h>

namespace Kpz
{
	struct DeviceSystemSize
	{
		int m_mDimXW, m_mDimYW;
		int m_lDimXW, m_lDimYW;
		int m_blocksX;
		int m_wLinesBlock, m_modBlockLines;
		int m_mThreads;

		__host__ std::ostream& print(std::ostream& o) {
			o << "[DBG] kpzSCA DeviceSystemSize:\n"
				"[DBG]\tSystem dimensions (quad): " << m_mDimXW+1 << " " << m_mDimYW+1
				<< "\n[DBG]\tblocksX = " << m_blocksX
				<< "\n[DBG]\twLinesBlock= " << m_wLinesBlock << " ; modBlockLines= " << m_modBlockLines << '\n';
			return o;
		}
	};

	__constant__ Kpz::DeviceSystemSize deviceSystemSize;
	/*!
	 * Also calculates grid/blocks dimx
	 */
	template<class GpuRng>
	void SchedulerSCA<GpuRng>::initStaticDeviceConstDyn(const Kpz::SystemSize& size)
	{
		DeviceSystemSize tmp;
		tmp.m_mDimXW = size.mDimXW();
		tmp.m_mDimYW = size.mDimYW();
		tmp.m_lDimXW = size.lDimXW();
		tmp.m_lDimYW = size.lDimYW();

		int lMaxThreadsArch = Kmc::log2((unsigned)m_prop->maxThreadsPerBlock);
		assert(lMaxThreadsArch <= Kpz::L_SCA_MAX_THREADS);

		//! \todo allow spanning block
		m_lThreads = std::min(tmp.m_lDimXW-1, lMaxThreadsArch);
		tmp.m_mThreads = (1<<m_lThreads)-1;
		m_blocksX = tmp.m_blocksX = 1<<(tmp.m_lDimXW-1-m_lThreads);
		m_blocks = m_prop->multiProcessorCount << (lMaxThreadsArch-m_lThreads);
		if((m_prop->major == 3 && (m_prop->minor == 7 | m_prop->minor == 5)) || m_prop->major == 5)
		{
			m_blocks <<= 1;
			std::cout << "sm_37,35 or 5* detected, doubling number of blocks." << std::endl;
		}
		// 1d checkerboard in y
		tmp.m_wLinesBlock = (1<<(tmp.m_lDimYW-1))/(m_blocks);
		if(tmp.m_wLinesBlock<4) // magic guess
		{
			m_blocks = m_prop->multiProcessorCount;
			tmp.m_wLinesBlock = (1<<(tmp.m_lDimYW-1))/(m_blocks);
			if(tmp.m_wLinesBlock < 1)
			{
				m_blocks = (1<<(tmp.m_lDimYW-1));
				tmp.m_wLinesBlock = 1;
			}
		}

		tmp.m_modBlockLines = (1<<(tmp.m_lDimYW-1))%(m_blocks);

		m_updatesPerGen = tmp.m_wLinesBlock;
		if(tmp.m_modBlockLines > 0)
			m_updatesPerGen += 1; // round up for mod blocks
		// for each line: 8 updates per word * words per line / number of threads 
		m_updatesPerGen <<= (3+tmp.m_lDimXW-m_lThreads);

		tmp.print(std::cerr << "[DBG] Setting SCA DeviceSystemSize:\n") << "[DBG] updatesPerGen= " << m_updatesPerGen  << '\n';
		std::cout << "Running " << m_blocks << " blocks and " << (1<<m_lThreads) << " threads each.\n";
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(deviceSystemSize, (void*)&tmp, sizeof(tmp)));

		delete m_gpuRng;
		m_gpuRng = new GpuRng(m_blocks << m_lThreads, &m_dsfmt);
	}
}
#define KPZ_DEVICE_SYSTEM_SIZE Kpz::deviceSystemSize

#define PARITY_BLOCK (parity&1)
#define PARITY_LATTICE ((parity>>1)&1)

__device__ size_t d_index(size_t d_x, size_t d_y)
{
	return (d_x&KPZ_DEVICE_SYSTEM_SIZE.m_mDimXW)
		+ ((d_y&KPZ_DEVICE_SYSTEM_SIZE.m_mDimYW)<<KPZ_DEVICE_SYSTEM_SIZE.m_lDimXW);
}

template<int L_THREADS, class KpzCore>
__global__ void kpzSCAKernelTemplate_span(unsigned int* d_System, int parity, KpzCore core)
{
#if 1 // shared and local alloc
	__shared__ unsigned int transfer[(1<<L_THREADS)];
	__shared__ unsigned int transfer2[(1<<L_THREADS)];
	unsigned int localSystem[4];
	typename KpzCore::Device deviceCore(core, PARITY_LATTICE);
#endif
	
	size_t d_x = 0;
	size_t d_y = ((PARITY_BLOCK) ?
			(KPZ_DEVICE_SYSTEM_SIZE.m_wLinesBlock + (blockIdx.x < KPZ_DEVICE_SYSTEM_SIZE.m_modBlockLines)) : 0)
		+ ((blockIdx.x*KPZ_DEVICE_SYSTEM_SIZE.m_wLinesBlock +
			( (blockIdx.x) >= KPZ_DEVICE_SYSTEM_SIZE.m_modBlockLines ? KPZ_DEVICE_SYSTEM_SIZE.m_modBlockLines : blockIdx.x))<<1);

	// load first line
	localSystem[0] = d_System[d_index(d_x + (threadIdx.x<<1), d_y)];
	localSystem[1] = d_System[d_index(d_x + (threadIdx.x<<1)+1, d_y)];

	// loop over lines in block
	const int Y = d_y + KPZ_DEVICE_SYSTEM_SIZE.m_wLinesBlock
		+ (blockIdx.x < KPZ_DEVICE_SYSTEM_SIZE.m_modBlockLines);
	for(; d_y < Y; ++d_y)
	{
		// load 'second' line
		localSystem[2] = d_System[d_index(d_x + (threadIdx.x<<1), d_y+1)];
		localSystem[3] = d_System[d_index(d_x + (threadIdx.x<<1)+1, d_y+1)];

		deviceCore.initLine(d_index(d_x + (threadIdx.x<<1), d_y));
		// update
		deviceCore.update(localSystem);

		// shift data left
		transfer[threadIdx.x] = localSystem[0];
		transfer2[threadIdx.x] = localSystem[2];
		__syncthreads();
		localSystem[0] = localSystem[1];
		localSystem[1] = transfer[(threadIdx.x+1)&KPZ_DEVICE_SYSTEM_SIZE.m_mThreads];
		localSystem[2] = localSystem[3];
		localSystem[3] = transfer2[(threadIdx.x+1)&KPZ_DEVICE_SYSTEM_SIZE.m_mThreads];
		/*localSystem[3] = (threadIdx.x+1>>L_THREADS) ? 0 : transfer[threadIdx.x+1];*/
	
		deviceCore.shiftLine();
		// update
		deviceCore.update(localSystem);

		// shift data back (right)
		transfer[(threadIdx.x+1)&KPZ_DEVICE_SYSTEM_SIZE.m_mThreads] = localSystem[1];
		transfer2[(threadIdx.x+1)&KPZ_DEVICE_SYSTEM_SIZE.m_mThreads] = localSystem[3];
		__syncthreads();
		localSystem[1] = localSystem[0];
		localSystem[0] = transfer[threadIdx.x];
		localSystem[3] = localSystem[2];
		localSystem[2] = transfer2[threadIdx.x];

		// write 'first' line
		d_System[d_index(d_x + (threadIdx.x<<1), d_y)] = localSystem[0];
		d_System[d_index(d_x + (threadIdx.x<<1)+1, d_y)] = localSystem[1];

		// shift data up
		localSystem[0] = localSystem[2];
		localSystem[1] = localSystem[3];
	}
	// write last 'second', now 'first', line
	d_System[d_index(d_x + (threadIdx.x<<1), d_y)] = localSystem[0];
	d_System[d_index(d_x + (threadIdx.x<<1)+1, d_y)] = localSystem[1];
}

template<int L_THREADS, class KpzCore>
__global__ void kpzSCAKernelTemplate_vertical(unsigned int* d_System, int parity, KpzCore core)
{
#if 1 // shared and local alloc
	__shared__ unsigned int transfer[(1<<L_THREADS)+1];
	__shared__ unsigned int transfer2[(1<<L_THREADS)+1];
	unsigned int localSystem[4];
	typename KpzCore::Device deviceCore(core, PARITY_LATTICE);
#endif

	const int Y = ((PARITY_BLOCK) ?
			(KPZ_DEVICE_SYSTEM_SIZE.m_wLinesBlock + (blockIdx.x < KPZ_DEVICE_SYSTEM_SIZE.m_modBlockLines)) : 0)
		+ ((blockIdx.x*KPZ_DEVICE_SYSTEM_SIZE.m_wLinesBlock +
			( (blockIdx.x) >= KPZ_DEVICE_SYSTEM_SIZE.m_modBlockLines ? KPZ_DEVICE_SYSTEM_SIZE.m_modBlockLines : blockIdx.x))<<1)
				  + KPZ_DEVICE_SYSTEM_SIZE.m_wLinesBlock
		+ (blockIdx.x < KPZ_DEVICE_SYSTEM_SIZE.m_modBlockLines);

	// loop over columns (x)
	for(size_t d_x = 0; d_x <= KPZ_DEVICE_SYSTEM_SIZE.m_mDimXW; d_x += (blockDim.x<<1))
	{

	size_t d_y = ((PARITY_BLOCK) ?
			(KPZ_DEVICE_SYSTEM_SIZE.m_wLinesBlock + (blockIdx.x < KPZ_DEVICE_SYSTEM_SIZE.m_modBlockLines)) : 0)
		+ ((blockIdx.x*KPZ_DEVICE_SYSTEM_SIZE.m_wLinesBlock +
			( (blockIdx.x) >= KPZ_DEVICE_SYSTEM_SIZE.m_modBlockLines ? KPZ_DEVICE_SYSTEM_SIZE.m_modBlockLines : blockIdx.x))<<1);

	// load first line
	localSystem[0] = d_System[d_index(d_x + (threadIdx.x<<1), d_y)];
	localSystem[1] = d_System[d_index(d_x + (threadIdx.x<<1)+1, d_y)];
	// loop over lines in block
	for(; d_y < Y; ++d_y)
	{
		// load 'second' line
		localSystem[2] = d_System[d_index(d_x + (threadIdx.x<<1), d_y+1)];
		localSystem[3] = d_System[d_index(d_x + (threadIdx.x<<1)+1, d_y+1)];
		// get end of first line, needed later
		if(threadIdx.x == blockDim.x-1) transfer[threadIdx.x+1] = d_System[d_index(d_x+(blockDim.x<<1), d_y)];

		deviceCore.initLine(d_index(d_x + (threadIdx.x<<1), d_y));
		// update
		deviceCore.update(localSystem);

		// shift data left
		transfer[threadIdx.x] = localSystem[0];
		transfer2[threadIdx.x] = localSystem[2];
		__syncthreads();
		localSystem[0] = localSystem[1];
		localSystem[1] = transfer[(threadIdx.x+1)];
		localSystem[2] = localSystem[3];
		localSystem[3] = transfer2[(threadIdx.x+1)];
		/*localSystem[3] = (threadIdx.x+1>>L_THREADS) ? 0 : transfer[threadIdx.x+1];*/
	
		deviceCore.shiftLine();
		// update
		deviceCore.update(localSystem);

		// shift data back (right)
		transfer[(threadIdx.x+1)] = localSystem[1];
		transfer2[(threadIdx.x+1)] = localSystem[3];
		__syncthreads();
		localSystem[1] = localSystem[0];
		localSystem[0] = transfer[threadIdx.x];
		localSystem[3] = localSystem[2];
		localSystem[2] = transfer2[threadIdx.x];
		if(threadIdx.x == blockDim.x-1) d_System[d_index(d_x+(blockDim.x<<1), d_y)] = transfer[threadIdx.x+1];

		// write 'first' line
		d_System[d_index(d_x + (threadIdx.x<<1), d_y)] = localSystem[0];
		d_System[d_index(d_x + (threadIdx.x<<1)+1, d_y)] = localSystem[1];

		// shift data up
		localSystem[0] = localSystem[2];
		localSystem[1] = localSystem[3];
	}
	// write last 'second', now 'first', line
	d_System[d_index(d_x + (threadIdx.x<<1), d_y)] = localSystem[0];
	d_System[d_index(d_x + (threadIdx.x<<1)+1, d_y)] = localSystem[1];

	}
}

// begin SchedulerSCA::callers
#ifndef KPZ_SWITCH_PRAND
#include "kpzRandom.cu"
#else
#include "PrandDevice.cu"
#endif
#include "kpzSCACores.cu"
#include "scheduler.h"

#if USE_DYNAMIC_RND_CACHE
#include "GPURndCache.cuh"
#endif

template<class GpuRng>
void Kpz::SchedulerSCA<GpuRng>::caller_mcsPQSyncRng(int parity)
{
#ifdef USE_CACHED_RND_OPT
	if((m_disorderP[0] == .5 || m_disorderP[0] == 0) && (m_disorderQ[0] == .5 || m_disorderQ[0] == 0))
	{
#if USE_DYNAMIC_RND_CACHE
		KpzSCACores::PQSyncRng<BinaryRndFloatCache_5<GpuRng> > core(*reinterpret_cast<BinaryRndFloatCache_5<GpuRng>* >(m_gpuRng), m_updatesPerGen, m_disorderP[0], m_disorderQ[0]);
#else
		KpzSCACores::PQSyncRngCacheing_5<GpuRng> core(*m_gpuRng, m_updatesPerGen, m_disorderP[0], m_disorderQ[0]);
#endif
		if(m_blocksX > 1)
			kpzSCAKernelTemplate_vertical<Kpz::L_SCA_MAX_THREADS> <<< m_blocks, 1<<m_lThreads >>>
				(d_System, parity, core);
		else
			kpzSCAKernelTemplate_span<Kpz::L_SCA_MAX_THREADS> <<< m_blocks, 1<<m_lThreads >>>
				(d_System, parity, core);
	}
	else
#endif
	{
		KpzSCACores::PQSyncRng<GpuRng> core(*m_gpuRng, m_updatesPerGen, m_disorderP[0], m_disorderQ[0]);
		if(m_blocksX > 1)
			kpzSCAKernelTemplate_vertical<Kpz::L_SCA_MAX_THREADS> <<< m_blocks, 1<<m_lThreads >>>
				(d_System, parity, core);
		else
			kpzSCAKernelTemplate_span<Kpz::L_SCA_MAX_THREADS> <<< m_blocks, 1<<m_lThreads >>>
				(d_System, parity, core);
	}
}

template<class GpuRng>
void Kpz::SchedulerSCA<GpuRng>::initDisorder2DeviceConst()
{
	float tmp[2];
	for(int a = 0; a < 2; ++a)
		tmp[a] = m_disorderP[a];
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(KpzSCACores::DEVICE_CONST_DISORDER2_P, (void*)tmp, 8));
	for(int a = 0; a < 2; ++a)
		tmp[a] = m_disorderQ[a];
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(KpzSCACores::DEVICE_CONST_DISORDER2_Q, (void*)tmp, 8));
}

template<class GpuRng>
void Kpz::SchedulerSCA<GpuRng>::caller_mcsDisorder2SyncRng(int parity)
{
	KpzSCACores::Disorder2SyncRng<GpuRng> core(*m_gpuRng, d_Disorder, m_updatesPerGen);
	if(m_blocksX > 1)
		kpzSCAKernelTemplate_vertical<Kpz::L_SCA_MAX_THREADS> <<< m_blocks, 1<<m_lThreads >>>
			(d_System, parity, core);
	else
		kpzSCAKernelTemplate_span<Kpz::L_SCA_MAX_THREADS> <<< m_blocks, 1<<m_lThreads >>>
			(d_System, parity, core);
}
#include "schedulerSCAInstances.h"
