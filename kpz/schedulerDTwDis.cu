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

#include "schedulerDTwDis.h"

#include "kpzConst.h"
#include "kpz.h"

#include <utility>

#include <cuda.h>

template <class GpuRng, class LocalLayout>
Kpz::SchedulerDTwDis<GpuRng, LocalLayout>::SchedulerDT(int lx, int ly, int device)
	: Base(lx, ly, device), d_Disorder(0)
{
}

template <class GpuRng, class LocalLayout>
Kpz::SchedulerDTwDis<GpuRng, LocalLayout>::SchedulerDT(const SystemSize& size, int device)
	: Base(size, device), d_Disorder(0)
{
}

template <class GpuRng, class LocalLayout>
Kpz::SchedulerDTwDis<GpuRng, LocalLayout>::SchedulerDT(SchedulerService &&other, int device)
	: Base(std::move(other), device)
{
}

template <class GpuRng, class LocalLayout>
Kpz::SchedulerDTwDis<GpuRng, LocalLayout>::~SchedulerDT()
{
	releaseDeviceDisorder();
}

template <class GpuRng, class LocalLayout>
void Kpz::SchedulerDTwDis<GpuRng, LocalLayout>::releaseDeviceDisorder()
{
	if(d_Disorder)
	{
		CUDA_SAFE_CALL(cudaFree(d_Disorder));
		d_Disorder = 0;
	}
}

template <class GpuRng, class LocalLayout>
bool Kpz::SchedulerDTwDis<GpuRng, LocalLayout>::setSize(int lx, int ly)
{
	if(!Base::setSize(lx,ly))
		return false;
	releaseDeviceDisorder();
	return true;
}

template <class GpuRng, class LocalLayout>
void Kpz::SchedulerDTwDis<GpuRng, LocalLayout>::mcsDisorder2SyncRng(int n)
{
	if(!d_Disorder)
	{
		CUDA_SAFE_CALL(cudaSetDevice(m_device));
		CUDA_SAFE_CALL(cudaMalloc((void**)&d_Disorder, m_size.sizeW()<<2));
		cudaError_t error = (cudaMemcpy((void*)d_Disorder, (void*)m_disorder, m_size.sizeW()<<2, cudaMemcpyHostToDevice)); 
		if(error != cudaSuccess)
			throw KmcExcept::CUDAError(error, __FILE__, __LINE__);
	}
	doMcs(n, &SchedulerDTwDis::caller_mcsDisorder2SyncRng);
}

#include "schedulerDTwDisInstances.h"
