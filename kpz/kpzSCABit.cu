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

#include "schedulerSCABit.h"

#include "kpzConst.h"

#include <kmcMath.h>

#include <iostream>
#include <algorithm>
#include <cassert>

#include <cudaError.h>

namespace Kpz
{
	struct DeviceSystemSizeBit
	{
		int m_lDimX, m_lDimY;
		size_t m_quadrant;
		int m_linesBlock, m_modBlockLines;

		__device__ inline int dimX() const {return 1<<m_lDimX;}
		__device__ inline int dimY() const {return 1<<m_lDimY;}

		__host__ std::ostream& print(std::ostream& o) {
			o << "[DBG] kpzSCABit DeviceSystemSize:\n"
				"[DBG]\tSystem dimensions (log2): " << m_lDimX << " " << m_lDimY
				<< "\n[DBG]\tquadrant = " << m_quadrant
				<< "\n[DBG]\tlinesBlock= " << m_linesBlock << " ; modBlockLines= " << m_modBlockLines << '\n';
			return o;
		}
	};

	__constant__ Kpz::DeviceSystemSizeBit deviceSystemSizeBit;

	const int L_VECTOR_SIZE_SITES = 5;
	const int L_WARP_SIZE = 5;
	const int M_WARP_SIZE = (1<<L_WARP_SIZE)-1;
	const int L_SCABIT_MAX_THREADS = 10;
	/*!
	 * Also calculates grid/blocks dimx
	 */
	template<class GpuRng>
	void SchedulerSCABit<GpuRng>::initStaticDeviceConstDyn(const Kpz::SystemSize& size)
	{
		assert(L_WARP_SIZE <= size.lDimX()-6);
		DeviceSystemSizeBit tmp;
		tmp.m_lDimX = size.lDimX();
		tmp.m_lDimY = size.lDimY();
		tmp.m_quadrant = size.sizeW()>>2;

		const int lMaxThreadsPerBlock = Kmc::log2((unsigned)m_prop->maxThreadsPerBlock);
		m_lThreads = std::min(std::min(L_SCABIT_MAX_THREADS, lMaxThreadsPerBlock), size.lDimX()-6);
		const int lBlockPerMPArch = std::min(lMaxThreadsPerBlock - m_lThreads, 16); // NVIDIA devices cannot scheduler more than 16 blocks per MP
		assert(lBlockPerMPArch >= 0);

		m_blocks = m_prop->multiProcessorCount << (lBlockPerMPArch);
		/*std::cerr << "sm_" << m_prop->major << m_prop->minor << '\n';*/
		if((m_prop->major == 3 && (m_prop->minor == 7 | m_prop->minor == 5)) || m_prop->major == 5)
		{
			m_blocks <<= 1;
			std::cout << "sm_37,35 or 5* detected, doubling number of blocks." << std::endl;
			if(m_prop->major == 5)
				std::cout << "... sm_50: doubled number of blocks not optimal if p=0.5. TODO: adjust dynamically." << std::endl;
		}

		tmp.m_linesBlock = (1<<(tmp.m_lDimY))/(m_blocks);
		if(tmp.m_linesBlock < 4)
		{
			m_blocks = ((1<<(tmp.m_lDimY))/m_prop->multiProcessorCount/4)*m_prop->multiProcessorCount;
			tmp.m_linesBlock = (1<<(tmp.m_lDimY))/(m_blocks);
		}
		tmp.m_modBlockLines = (1<<(tmp.m_lDimY))%(m_blocks);

		m_updatesPerGen = tmp.m_linesBlock;
		if(tmp.m_modBlockLines > 0)
			m_updatesPerGen += 1; // round up for mod blocks
		m_updatesPerGen <<= tmp.m_lDimX-m_lThreads;

		std::cout << "Running " << m_blocks << " blocks and " << (1<<m_lThreads) << " threads each.\n";
		tmp.print(std::cerr);
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(deviceSystemSizeBit, (void*)&tmp, sizeof(tmp)));

		delete m_gpuRng;
		m_gpuRng = new GpuRng(m_blocks << m_lThreads, &m_dsfmt);
	}
#define KPZ_DEVICE_SYSTEM_SIZE Kpz::deviceSystemSizeBit

#define PARITY_LATTICE (parity)
#define PARITY_LINE ( PARITY_LATTICE == (d_y&1) )

template<class KpzCore>
__global__ void kpzSCABitKernelTemplate_span(unsigned int* d_System, int parity, KpzCore core)
{
#if 1 // shared and local alloc
	typename KpzCore::Device deviceCore(core);
	__shared__ unsigned int sxDataX[1<<(L_SCABIT_MAX_THREADS+1)];
#endif

	int d_y =
		((blockIdx.x*KPZ_DEVICE_SYSTEM_SIZE.m_linesBlock +
			( (blockIdx.x >= KPZ_DEVICE_SYSTEM_SIZE.m_modBlockLines) ? KPZ_DEVICE_SYSTEM_SIZE.m_modBlockLines : blockIdx.x)));

	// loop over lines in block
	const int Y = d_y + KPZ_DEVICE_SYSTEM_SIZE.m_linesBlock
		+ (blockIdx.x < KPZ_DEVICE_SYSTEM_SIZE.m_modBlockLines);
	for(; d_y < Y; ++d_y)
	{
		const size_t offset = (d_y << (KPZ_DEVICE_SYSTEM_SIZE.m_lDimX-Kpz::L_VECTOR_SIZE_SITES-1)) + threadIdx.x;
		size_t index = ((KPZ_DEVICE_SYSTEM_SIZE.m_quadrant<<1)&(0-PARITY_LATTICE))+offset;
		unsigned int thisDataX = d_System[index];
		unsigned int thisDataY = d_System[index+KPZ_DEVICE_SYSTEM_SIZE.m_quadrant];
		index = ((KPZ_DEVICE_SYSTEM_SIZE.m_quadrant<<1)&(~(0-PARITY_LATTICE)));
		const size_t yIndexY = index + KPZ_DEVICE_SYSTEM_SIZE.m_quadrant
			+ threadIdx.x+ (((d_y+1)&(KPZ_DEVICE_SYSTEM_SIZE.dimY()-1))<<(KPZ_DEVICE_SYSTEM_SIZE.m_lDimX-(Kpz::L_VECTOR_SIZE_SITES+1)));
		unsigned int yDataY = d_System[yIndexY];
		index += offset;
		sxDataX[threadIdx.x] = d_System[index];
		unsigned int xDataX = sxDataX[threadIdx.x];

		if(!PARITY_LINE)
		{
			xDataX >>= 1;
			__syncthreads();
			xDataX |= sxDataX[(threadIdx.x+1)&(blockDim.x-1)]<<31;
		}

		{
			typename KpzCore::Device::UpdateMask update(deviceCore);
			update.genMask(thisDataX, thisDataY, xDataX, yDataY);

			thisDataX ^= update.mask();
			thisDataY ^= update.mask();
			yDataY ^= update.mask();
			if(!PARITY_LINE)
			{
				sxDataX[threadIdx.x] ^= (update.mask()<<1);
				__syncthreads();
				sxDataX[(threadIdx.x+1)&(blockDim.x-1)] ^= (update.mask()>>31);
				/*atomicXor(&sxDataX[threadIdx.x], (update.mask()<<1));*/
				/*atomicXor(&sxDataX[(threadIdx.x+1)&(blockDim.x-1)], (update.mask()>>31));*/
				__syncthreads();
			}
			else
				sxDataX[threadIdx.x] ^= update.mask();
		}

		d_System[yIndexY] = yDataY;
		d_System[index] = sxDataX[threadIdx.x];
		index = ((KPZ_DEVICE_SYSTEM_SIZE.m_quadrant<<1)&(0-PARITY_LATTICE))+offset;
		d_System[index] = thisDataX;
		d_System[index+KPZ_DEVICE_SYSTEM_SIZE.m_quadrant] = thisDataY;
		/*printf("tidx %d :index= %d ,yIndexY= %d ,thisDataX= %X ,thisDataY= %X ,xDataXnext= %X ,yDataY= %X\n", threadIdx.x, (int)index, (int)yIndexY, thisDataX, thisDataY, xDataXnext, yDataY);*/
	}
}

const unsigned int M_WRAP_BUFFER_SIZE = (1<<(L_SCABIT_MAX_THREADS+1))-1;

template<class KpzCore>
__global__ void kpzSCABitKernelTemplate(unsigned int* d_System, int parity, KpzCore core)
{
#if 1 // shared and local alloc
	typename KpzCore::Device deviceCore(core);
	__shared__ unsigned int sxDataX[M_WRAP_BUFFER_SIZE+1];
#endif

	int d_y =
		((blockIdx.x*KPZ_DEVICE_SYSTEM_SIZE.m_linesBlock +
			( (blockIdx.x >= KPZ_DEVICE_SYSTEM_SIZE.m_modBlockLines) ? KPZ_DEVICE_SYSTEM_SIZE.m_modBlockLines : blockIdx.x)));

	// loop over lines in block
	const int Y = d_y + KPZ_DEVICE_SYSTEM_SIZE.m_linesBlock
		+ (blockIdx.x < KPZ_DEVICE_SYSTEM_SIZE.m_modBlockLines);
	for(; d_y < Y; ++d_y)
	{
		if(PARITY_LINE)
		{
			for(int d_xw = threadIdx.x; d_xw < (KPZ_DEVICE_SYSTEM_SIZE.dimX()>>6); d_xw += blockDim.x)
			{
				const size_t offset = (d_y << (KPZ_DEVICE_SYSTEM_SIZE.m_lDimX-L_VECTOR_SIZE_SITES-1)) + d_xw;
				size_t index = ((KPZ_DEVICE_SYSTEM_SIZE.m_quadrant<<1)&(0-PARITY_LATTICE))+offset;
				unsigned int thisDataX = d_System[index];
				unsigned int thisDataY = d_System[index+KPZ_DEVICE_SYSTEM_SIZE.m_quadrant];
				index = ((KPZ_DEVICE_SYSTEM_SIZE.m_quadrant<<1)&(~(0-PARITY_LATTICE)));
				const size_t yIndexY = index + KPZ_DEVICE_SYSTEM_SIZE.m_quadrant
					+ d_xw+ (((d_y+1)&(KPZ_DEVICE_SYSTEM_SIZE.dimY()-1))<<(KPZ_DEVICE_SYSTEM_SIZE.m_lDimX-(Kpz::L_VECTOR_SIZE_SITES+1)));
				unsigned int yDataY = d_System[yIndexY];
				index += offset;
				unsigned int xDataX = d_System[index];
				{
					typename KpzCore::Device::UpdateMask update(deviceCore);
					update.genMask(thisDataX, thisDataY, xDataX, yDataY);

					thisDataX ^= update.mask();
					thisDataY ^= update.mask();
					yDataY ^= update.mask();
					xDataX ^= update.mask();
				}

				d_System[yIndexY] = yDataY;
				d_System[index] = xDataX;
				index = ((KPZ_DEVICE_SYSTEM_SIZE.m_quadrant<<1)&(0-PARITY_LATTICE))+offset;
				d_System[index] = thisDataX;
				d_System[index+KPZ_DEVICE_SYSTEM_SIZE.m_quadrant] = thisDataY;
			}
		}
		else
		{
			sxDataX[threadIdx.x] = d_System[((KPZ_DEVICE_SYSTEM_SIZE.m_quadrant<<1)&(~(0-PARITY_LATTICE)))
				+ (d_y << (KPZ_DEVICE_SYSTEM_SIZE.m_lDimX-Kpz::L_VECTOR_SIZE_SITES-1)) + threadIdx.x];
			__shared__ unsigned int sxD0Cache;
			for(int d_xw = threadIdx.x; d_xw < (KPZ_DEVICE_SYSTEM_SIZE.dimX()>>6); )
			{
				const size_t offset = (d_y << (KPZ_DEVICE_SYSTEM_SIZE.m_lDimX-L_VECTOR_SIZE_SITES-1)) + d_xw;
				size_t index = ((KPZ_DEVICE_SYSTEM_SIZE.m_quadrant<<1)&(0-PARITY_LATTICE))+offset;
				unsigned int thisDataX = d_System[index];
				unsigned int thisDataY = d_System[index+KPZ_DEVICE_SYSTEM_SIZE.m_quadrant];
				index = ((KPZ_DEVICE_SYSTEM_SIZE.m_quadrant<<1)&(~(0-PARITY_LATTICE)));
				const size_t yIndexY = index + KPZ_DEVICE_SYSTEM_SIZE.m_quadrant
					+ d_xw+ (((d_y+1)&(KPZ_DEVICE_SYSTEM_SIZE.dimY()-1))<<(KPZ_DEVICE_SYSTEM_SIZE.m_lDimX-(Kpz::L_VECTOR_SIZE_SITES+1)));
				unsigned int yDataY = d_System[yIndexY];
				index += offset + blockDim.x;
				/*unsigned int xDataX = sxDataX[d_xw&M_WRAP_BUFFER_SIZE];*/

				if(d_xw + blockDim.x < (KPZ_DEVICE_SYSTEM_SIZE.dimX()>>6))
				{
					sxDataX[(d_xw+blockDim.x)&M_WRAP_BUFFER_SIZE] = d_System[index];
				}
				else if(threadIdx.x == 0)
					sxDataX[(d_xw+blockDim.x)&M_WRAP_BUFFER_SIZE] = sxD0Cache; // place cached first word after end


				/*xDataX >>= 1;*/
				/*__syncthreads();*/
				/*xDataX |= sxDataX[(d_xw+1)&M_WRAP_BUFFER_SIZE]<<31;*/

				{
					// WARNING: inline mask calc was replaced by KpzCore::Device::UpdateMask here, but was not tested
					typename KpzCore::Device::UpdateMask update(deviceCore);
					__syncthreads();
					// random was before sync, mask after
					update.genMask(thisDataX, thisDataY, (sxDataX[d_xw&M_WRAP_BUFFER_SIZE]>>1) | (sxDataX[(d_xw+1)&M_WRAP_BUFFER_SIZE]<<31), yDataY);

					thisDataX ^= update.mask();
					thisDataY ^= update.mask();
					yDataY ^= update.mask();
					__syncthreads();
					sxDataX[d_xw&M_WRAP_BUFFER_SIZE] ^= update.mask()<<1;
					__syncthreads();
					sxDataX[(d_xw+1)&M_WRAP_BUFFER_SIZE] ^= update.mask()>>31;
				}

				d_System[yIndexY] = yDataY;
				__syncthreads();
				if(d_xw > 0)
					d_System[index - blockDim.x] = sxDataX[d_xw&M_WRAP_BUFFER_SIZE];
				else
					sxD0Cache = sxDataX[0];
				index = ((KPZ_DEVICE_SYSTEM_SIZE.m_quadrant<<1)&(0-PARITY_LATTICE))+offset;
				d_System[index] = thisDataX;
				d_System[index+KPZ_DEVICE_SYSTEM_SIZE.m_quadrant] = thisDataY;
				d_xw += blockDim.x;
				if(d_xw  == (KPZ_DEVICE_SYSTEM_SIZE.dimX()>>6)) // threadIdx.x == 0
					d_System[((KPZ_DEVICE_SYSTEM_SIZE.m_quadrant<<1)&(~(0-PARITY_LATTICE)))
						+ (d_y << (KPZ_DEVICE_SYSTEM_SIZE.m_lDimX-Kpz::L_VECTOR_SIZE_SITES-1))]
						= sxDataX[(d_xw)&M_WRAP_BUFFER_SIZE]; // store first word
			}
		}
	}
}

} // namespace Kpz

#ifndef KPZ_SWITCH_PRAND
#include "kpzRandom.cu"
#else
#include "PrandDevice.cu"
#endif
#include "kpzSCABitCores.cu"
#include "scheduler.h"

template<class GpuRng>
template<class KpzCore>
void Kpz::SchedulerSCABit<GpuRng>::run(int parity, KpzCore core)
{
	if(m_size.dimX() == (1<<m_lThreads<<6))
	{
		kpzSCABitKernelTemplate_span <<< m_blocks, (1<<m_lThreads) >>>
			(d_System, parity, core);
	}
	else
	{
		kpzSCABitKernelTemplate <<< m_blocks, (1<<m_lThreads) >>>
			(d_System, parity, core);
	}
}

template<class GpuRng>
void Kpz::SchedulerSCABit<GpuRng>::caller_mcsSyncRng(int parity)
{
	run(parity, KpzSCABitCores::Empty<GpuRng>(*m_gpuRng, m_updatesPerGen));
}

template<class GpuRng>
__device__ inline unsigned int genP0_875(KpzSCABitCores::Base<typename GpuRng::Device>& dev)
{
	return dev.random() | dev.random() | dev.random();
}

template<class GpuRng>
__device__ inline unsigned int genP0_75(KpzSCABitCores::Base<typename GpuRng::Device>& dev)
{
	return dev.random() | dev.random();
}

template<class GpuRng>
__device__ inline unsigned int genP0_625(KpzSCABitCores::Base<typename GpuRng::Device>& dev)
{
	return (dev.random() & dev.random()) | dev.random();
}

template<class GpuRng>
__device__ inline unsigned int genP0_375(KpzSCABitCores::Base<typename GpuRng::Device>& dev)
{
	return (dev.random() | dev.random()) & dev.random();
}

template<class GpuRng>
__device__ inline unsigned int genP0_25(KpzSCABitCores::Base<typename GpuRng::Device>& dev)
{
	return dev.random() & dev.random();
}

template<class GpuRng>
__device__ inline unsigned int genP0_125(KpzSCABitCores::Base<typename GpuRng::Device>& dev)
{
	return dev.random() & dev.random() & dev.random();
}

template<class GpuRng>
__device__ inline unsigned int genP1_32(KpzSCABitCores::Base<typename GpuRng::Device>& dev)
{
	return dev.random() & dev.random() & dev.random() & dev.random() & dev.random();
}

template<class GpuRng>
void Kpz::SchedulerSCABit<GpuRng>::caller_mcsPQSyncRng(int parity)
{
	if(m_disorderQ[0]==0.)
	{
		if(m_disorderP[0] == .75)
			run(parity, KpzSCABitCores::EmptyGenP<GpuRng, genP0_75<GpuRng> >(*m_gpuRng, m_updatesPerGen));
		else
			run(parity, KpzSCABitCores::EmptyAnyP<GpuRng>(*m_gpuRng, m_updatesPerGen, m_disorderP[0]));
	}
	else if(m_disorderQ[0] == m_disorderP[0])
	{
		if(m_disorderQ[0] == .5)
		{
			run(parity, KpzSCABitCores::EmptyPQ<GpuRng>(*m_gpuRng, m_updatesPerGen));
		}
		else if(m_disorderQ[0] == 0.03125)
		{
			run(parity, KpzSCABitCores::EmptyGenPQ<GpuRng, genP1_32<GpuRng> >(*m_gpuRng, m_updatesPerGen));
		}
		else if(m_disorderQ[0] == 0.75)
		{
			run(parity, KpzSCABitCores::EmptyGenPQ<GpuRng, genP0_75<GpuRng> >(*m_gpuRng, m_updatesPerGen));
		}
		else if(m_disorderQ[0] == 0.875)
		{
			run(parity, KpzSCABitCores::EmptyGenPQ<GpuRng, genP0_875<GpuRng> >(*m_gpuRng, m_updatesPerGen));
		}
		else if(m_disorderQ[0] == 0.625)
		{
			run(parity, KpzSCABitCores::EmptyGenPQ<GpuRng, genP0_625<GpuRng> >(*m_gpuRng, m_updatesPerGen));
		}
		else if(m_disorderQ[0] == 0.375)
		{
			run(parity, KpzSCABitCores::EmptyGenPQ<GpuRng, genP0_375<GpuRng> >(*m_gpuRng, m_updatesPerGen));
		}
		else if(m_disorderQ[0] == 0.25)
		{
			run(parity, KpzSCABitCores::EmptyGenPQ<GpuRng, genP0_25<GpuRng> >(*m_gpuRng, m_updatesPerGen));
		}
		else if(m_disorderQ[0] == 0.125)
		{
			run(parity, KpzSCABitCores::EmptyGenPQ<GpuRng, genP0_125<GpuRng> >(*m_gpuRng, m_updatesPerGen));
		}
	}
	else
	{
		run(parity, KpzSCABitCores::EmptyAnyPQ<GpuRng>(*m_gpuRng, m_updatesPerGen, m_disorderP[0], m_disorderQ[0]));
	}
}

#include "schedulerSCABitInstances.h"
