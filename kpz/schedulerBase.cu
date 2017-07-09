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

#include "schedulerBase.h"

#include "kpzConst.h"
#include "kpz.h"
#include "systemSize.h"

#include <kmcExceptCUDA.h>
#include <benchmarking.h>

#include <cuda.h>
#include <cudaError.h>
#include <iomanip>

void Kpz::SchedulerBase::pushSystem()
{
	join();
	cudaError_t error = (cudaMemcpy((void*)d_System, (void*)m_system, m_size.sizeW()<<2, cudaMemcpyHostToDevice)); 
	if(error == cudaSuccess)
		m_deviceMcs = m_mcs;
	else
		throw KmcExcept::CUDAError(error, __FILE__, __LINE__);
}

void Kpz::SchedulerBase::popSystem()
{
	join();
	cudaError_t error = (cudaMemcpy((void*)m_system, (void*)d_System, m_size.sizeW()<<2, cudaMemcpyDeviceToHost)); 
	if(error == cudaSuccess)
		m_mcs = m_deviceMcs;
	else
		throw KmcExcept::CUDAError(error, __FILE__, __LINE__);
	#if 0
	std::cerr << std::hex << '\n' << m_system[0] << ' ' << m_system[1] 
		//<< '\n' << m_system[2] << ' ' << m_system[3]
		<< '\n' << m_system[4] << ' ' << m_system[5] << '\n';
	std::cerr << '\n';
	for(int a = 0; a < THREADS; ++a)
	{
		std::cerr << "Thread " << a << ": " << m_system[4*a] << ", " << m_system[4*a+1] << ", " << m_system[4*a+2] << '\n';
	}
	#endif
}

void Kpz::SchedulerBase::initCUDA()
{
	CUDA_SAFE_CALL(cudaMalloc((void**)&d_System, m_size.sizeW()<<2));
	if(!d_System)
	{
		std::cerr << "Insufficient devicememory.\n";
		exit(1);
	}
	CUDA_SAFE_CALL(cudaHostRegister(m_system, m_size.sizeW()<<2, cudaHostRegisterPortable));
}

Kpz::SchedulerBase::~SchedulerBase()
{
	if (d_System)
		CUDA_SAFE_CALL(cudaFree(d_System));
	CUDA_SAFE_CALL(cudaHostUnregister(m_system));
}

void Kpz::SchedulerBase::setMcs(unsigned int mcs)
{
	SchedulerService::setMcs(mcs);
	m_deviceMcs = m_mcs;
}

// async mcs* calls
void* Kpz::SchedulerBase::callMcs(void* vargs)
{
	timeAinit;
	timeAstart;
	auto args = reinterpret_cast<CallMcsArgs*>(vargs);
	CUDA_SAFE_CALL(cudaSetDevice(dynamic_cast<SchedulerBase*>(args->scheduler)->m_device));
	(args->scheduler->*(args->fkt))(args->n);
	timeAstopS("SchedulerBase::callMcs_" << args->n);
	return 0;
}
