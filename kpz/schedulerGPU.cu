/***************************************************************************
*   Copyright 2011 - 2012 Jeffrey Kelling <j.kelling@hzdr.de>
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

#include "schedulerGPU.h"

#include <cuda.h>
#include <cudaError.h>
#include <iostream>

int getSetDevice(int begin = 0)
{
	int ndev = 0;
	CUDA_SAFE_CALL(cudaGetDeviceCount(&ndev));
	cudaError_t err;
	for(int a = begin; a < ndev; ++a)
	{
		 /* cudaSetDevice returns cudaSuccess even when the device is in use in
		  * exclusive_process mode, ... */
		err = cudaSetDevice(a);
		/* ... cudaThreadSynchronize reliably fails when device is in use by
		 * another proces in exclusive_process mode. */
		err = cudaThreadSynchronize();
		if(err == cudaSuccess)
			return a;
		cudaGetLastError();
	}
	std::cerr << "CUDA ERROR at " << __FILE__ << ':' << __LINE__ << ": " << cudaGetErrorString(err)
		<< "[EE][SchedulerGPU] could not set device.\n";
	exit(1);
}

Kpz::SchedulerGPU::SchedulerGPU (int device, int smemWperBlock)
	: m_prop(new cudaDeviceProp), m_device(device)
{
	if(device < 0)
		m_device = getSetDevice();
	else
		m_device = getSetDevice(device);
	CUDA_SAFE_CALL(cudaGetDeviceProperties(m_prop, m_device));
	m_blocks = m_prop->multiProcessorCount;
	std::cout << "Using " << m_blocks << " Multiprocessors on device " << m_device << '.' << std::endl;
	if(smemWperBlock > 0)
		m_blocks = m_blocks * ((m_prop->sharedMemPerBlock>>2) / smemWperBlock);
	if(m_prop->major == 3 && m_prop->minor == 7)
	{
		m_blocks <<= 1;
		std::cout << "sm_37 detected, doubling number of blocks." << std::endl;
	}
	std::cout << "Running " << m_blocks << " blocks on device " << m_device << '.' << std::endl;
}

Kpz::SchedulerGPU::~SchedulerGPU ()
{
	delete m_prop;
}
