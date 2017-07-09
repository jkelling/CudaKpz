/***************************************************************************
*   Copyright 2014 - 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#include "schedulerDT.h"

#include "kpzConst.h"
#include "kpz.h"
#include "systemSize.h"
#include "LocalLayoutIO.h"

#include <shuffle.h>
#include <timer.h>
#include <kmcExceptCUDA.h>

#include <utility>

#include <cuda.h>
#include <cudaError.h>

template <class GpuRng, class LocalLayout>
void Kpz::SchedulerDT<GpuRng, LocalLayout>::doMcs(int n, Kpz::SchedulerDT<GpuRng, LocalLayout>::kernelCallerFct mcs)
{
	cudaEvent_t event;
	cudaEventCreate(&event);

	if(n <= 0)
		return;
	doMcs(mcs);
	for(--n; n > 0; --n)
	{
		cudaEventRecord(event);
		doMcs(mcs);
		while(cudaEventQuery(event) == cudaErrorNotReady);
		if(Kmc::TimerSingleton::atEnd())
			break;
	}
	cudaThreadSynchronize();
	cudaEventDestroy(event);
	cudaError_t error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		throw KmcExcept::CUDAError(error, __FILE__, __LINE__);
	}
}

template <class GpuRng, class LocalLayout>
void Kpz::SchedulerDT<GpuRng, LocalLayout>::doMcs(Kpz::SchedulerDT<GpuRng, LocalLayout>::kernelCallerFct mcs)
{
	for(int a = 0; a < (1<<MyLocalLayout::L_MCS_DIV); ++a)
	{
		const int xw = (int)(dsfmt_genrand_close_open(&m_dsfmt)*(2<<MyLocalLayout::L_BLOCK_DIM_X));
		const int yw = (int)(dsfmt_genrand_close_open(&m_dsfmt)*(2<<MyLocalLayout::L_BLOCK_DIM_Y));
		Kmc::Shuffle<4> seq(&m_dsfmt);
		for(int a = 0; a < 4; ++a)
		{
			(*mcs)(this, xw,yw, seq[a]);
		}
	}
	++m_deviceMcs;
}

template <class GpuRng, class LocalLayout>
void Kpz::SchedulerDT<GpuRng, LocalLayout>::init()
{
	std::cout << "SchedulerDT (GPU, DT ro)\n";
	std::cout << "\tBlock layer DD: " << INNER_DD_STRING << " -- " << INNER_DD_STRING_LONG << '\n';
	MyLocalLayout::print(std::cout);
	m_size.setGPUfromLocalLayout<MyLocalLayout>();
	m_size.adjustDT();
	initStaticDeviceConstDyn(m_size());
	m_randomNumberBlockMultiplier = calcRandomNumberBlockMultiplier();
	std::cout << "\nSet system size to " << m_size.lDimX() << ',' << m_size.lDimY() << '\n' << m_size;
}

template <class GpuRng, class LocalLayout>
Kpz::SchedulerDT<GpuRng, LocalLayout>::SchedulerDT(int lx, int ly, int device)
	: SchedulerBase(lx, ly, device, MyLocalLayout::BLOCK_DATA)
	, m_gpuRng(getRandomNumberGeneratorCount(MyLocalLayout::L_THREADS), &m_dsfmt)
	, d_Disorder(0)
{
	init();
}

template <class GpuRng, class LocalLayout>
Kpz::SchedulerDT<GpuRng, LocalLayout>::SchedulerDT(const SystemSize& size, int device)
	: SchedulerBase(size, device, MyLocalLayout::BLOCK_DATA)
	, m_gpuRng(getRandomNumberGeneratorCount(MyLocalLayout::L_THREADS), &m_dsfmt)
	, d_Disorder(0)
{
	init();
}

template <class GpuRng, class LocalLayout>
Kpz::SchedulerDT<GpuRng, LocalLayout>::SchedulerDT(SchedulerService &&other, int device)
	: SchedulerBase(std::move(other), device, MyLocalLayout::BLOCK_DATA)
	, m_gpuRng(getRandomNumberGeneratorCount(MyLocalLayout::L_THREADS), &m_dsfmt)
	, d_Disorder(0)
{
	init();
}

template <class GpuRng, class LocalLayout>
Kpz::SchedulerDT<GpuRng, LocalLayout>::~SchedulerDT()
{
	releaseDeviceDisorder();
}

template <class GpuRng, class LocalLayout>
bool Kpz::SchedulerDT<GpuRng, LocalLayout>::changedSizeEvent()
{
	if(!SchedulerBase::changedSizeEvent())
		return false;
	m_size.setGPUfromLocalLayout<MyLocalLayout>();
	m_size.adjustDT();
	m_randomNumberBlockMultiplier = calcRandomNumberBlockMultiplier();
	initStaticDeviceConstDyn(m_size());
	releaseDeviceDisorder();
	return true;
}

template <class GpuRng, class LocalLayout>
void Kpz::SchedulerDT<GpuRng, LocalLayout>::releaseDeviceDisorder()
{
	if(d_Disorder)
	{
		CUDA_SAFE_CALL(cudaFree(d_Disorder));
		d_Disorder = 0;
	}
}

#include <KMCsplash.h>

template <class GpuRng, class LocalLayout>
bool Kpz::SchedulerDT<GpuRng, LocalLayout>::collectH5(splash::DataCollector* data, const std::string& prefix)
{
	if(!SchedulerBase::collectH5(data, prefix))
		return false;
	m_gpuRng.writeH5(data, m_mcs, prefix);
	return true;
}

template <class GpuRng, class LocalLayout>
bool Kpz::SchedulerDT<GpuRng, LocalLayout>::retrieveH5(splash::DataCollector* data, const std::string& prefix)
{
	if(!SchedulerBase::retrieveH5(data, prefix))
		return false;
	if(!m_gpuRng.readH5(data, m_mcs, prefix))
	{
		m_gpuRng.randomize();
		std::cout << "[SchedulerDT][retrieveH5][WW] Failed to restore GPU RNG state from file, randomizing. Check if dSFMT was restored.\n";
	}
	return true;
}

namespace Kpz
{

template<class GpuRng, class LocalLayout>
struct SchedulerDT_mcsDisorder2SyncRng
{
	static void f(SchedulerDT<GpuRng, LocalLayout>* pthis, int n)
	{
		std::cerr << "[BUG] mcsDisorder2SyncRng() not implemented.\n";
		exit(1);
	}
};

template<class GpuRng >
struct SchedulerDT_mcsDisorder2SyncRng<GpuRng, SchedulerDT_LocalLayoutDis>
{
	static void f(SchedulerDT<GpuRng, SchedulerDT_LocalLayoutDis>* pthis, int n)
	{
		pthis->doMcs(n, &SchedulerDTCallers::caller_mcsDisorder2SyncRng<SchedulerDT<GpuRng, SchedulerDT_LocalLayoutDis> >);
	}
};

}

template <class GpuRng, class LocalLayout>
void Kpz::SchedulerDT<GpuRng, LocalLayout>::mcsDisorder2SyncRng(int n)
{
	if(!d_Disorder)
	{
		std::cerr << "alloc d_Disorder\n";
		CUDA_SAFE_CALL(cudaSetDevice(m_device));
		CUDA_SAFE_CALL(cudaMalloc((void**)&d_Disorder, m_size.sizeW()<<2));
		if(!d_Disorder)
		{
			std::cerr << "Insufficient devicememory to store disorder2.\n";
			exit(1);
		}
		cudaError_t error = (cudaMemcpy((void*)d_Disorder, (void*)m_disorder, m_size.sizeW()<<2, cudaMemcpyHostToDevice)); 
		if(error != cudaSuccess)
			throw KmcExcept::CUDAError(error, __FILE__, __LINE__);
	}
	SchedulerDT_mcsDisorder2SyncRng<GpuRng, LocalLayout>::f(this, n);
}

#include "schedulerDTInstances.h"
