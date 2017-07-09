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
#include "kpz.h"
#include "systemSize.h"

#include <timer.h>
#include <kmcExceptCUDA.h>

#include <utility>

#include <cuda.h>
#include <cudaError.h>

template<class GpuRng>
double Kpz::SchedulerSCA<GpuRng>::mcsScale() const
{
	return m_disorderP[0]+m_disorderQ[0];
}

template <class GpuRng>
void Kpz::SchedulerSCA<GpuRng>::doMcs(int n, Kpz::SchedulerSCA<GpuRng>::kernelCallerFct mcs)
{
	cudaEvent_t event;
	cudaEventCreate(&event);

	if(n <= 0)
		return;
	for(int a = 0; a < 4; ++a)
		(this->*mcs)(a);
	++m_deviceMcs;
	for(--n; n > 0; --n)
	{
		cudaEventRecord(event);
		for(int a = 0; a < 4; ++a)
			(this->*mcs)(a);
		++m_deviceMcs;
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

template <class GpuRng>
void Kpz::SchedulerSCA<GpuRng>::init()
{
	std::cout << "SchedulerSCA (GPU)\n";
	initStaticDeviceConstDyn(m_size);
	setPQDefault();
}

template <class GpuRng>
void Kpz::SchedulerSCA<GpuRng>::releaseDeviceDisorder()
{
	if(d_Disorder)
	{
		CUDA_SAFE_CALL(cudaFree(d_Disorder));
		d_Disorder = 0;
	}
}

template <class GpuRng>
Kpz::SchedulerSCA<GpuRng>::SchedulerSCA(int lx, int ly, int device)
	: SchedulerBase(lx, ly, device), m_gpuRng(0), d_Disorder(0)
{
	init();
}

template <class GpuRng>
Kpz::SchedulerSCA<GpuRng>::SchedulerSCA(const SystemSize& size, int device)
	: SchedulerBase(size, device), m_gpuRng(0), d_Disorder(0)
{
	init();
}

template <class GpuRng>
Kpz::SchedulerSCA<GpuRng>::SchedulerSCA(SchedulerService &&other, int device)
	: SchedulerBase(std::move(other), device), m_gpuRng(0), d_Disorder(0)
{
	init();
}

template <class GpuRng>
Kpz::SchedulerSCA<GpuRng>::~SchedulerSCA()
{
	delete m_gpuRng;
	releaseDeviceDisorder();
}

template <class GpuRng>
bool Kpz::SchedulerSCA<GpuRng>::changedSizeEvent()
{
	if(!SchedulerBase::changedSizeEvent())
		return false;
	initStaticDeviceConstDyn(m_size);
	releaseDeviceDisorder();
	return true;
}

template <class GpuRng>
void Kpz::SchedulerSCA<GpuRng>::mcsDisorder2SyncRng(int n)
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
		initDisorder2DeviceConst();
	}
	doMcs(n, &SchedulerSCA::caller_mcsDisorder2SyncRng);
}

#include <KMCsplash.h>

template <class GpuRng>
bool Kpz::SchedulerSCA<GpuRng>::collectH5(splash::DataCollector* data, const std::string& prefix)
{
	if(!SchedulerBase::collectH5(data, prefix))
		return false;
	m_gpuRng->writeH5(data, m_mcs, prefix);
	return true;
}

template <class GpuRng>
bool Kpz::SchedulerSCA<GpuRng>::retrieveH5(splash::DataCollector* data, const std::string& prefix)
{
	if(!SchedulerBase::retrieveH5(data, prefix))
		return false;
	if(!m_gpuRng->readH5(data, m_mcs, prefix))
	{
		m_gpuRng->randomize();
		std::cout << "[SchedulerSCA][retrieveH5][WW] Failed to restore GPU RNG state from file, randomizing. Check if dSFMT was restored.\n";
	}
	return true;
}

#include "schedulerSCAInstances.h"
