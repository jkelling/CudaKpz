/***************************************************************************
*   Copyright 2011 - 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#include "kpzFrontend.h"

#include <iostream>
#include <sstream>

#include <cuda.h>
#include <cudaError.h>
#include <mpi.h>

int Kpz::SchedulerConfig::initCUDA() const
{
	int device = -1;
	if(m_setGpuByMpiRank >= 0)
	{
		if(!MPI::Is_initialized())
		{
			std::cerr << "[WW][SchedulerConfig][initCUDA] MPI was not initialized before creating GPU scheduler with MPI flags."
				"Calling MPI::Init()\n";
			MPI::Init();
		}
		unsigned int rank = MPI::COMM_WORLD.Get_rank();
		device = m_setGpuByMpiRank+rank;
		CUDA_SAFE_CALL( cudaSetDevice(device) );
		std::ostringstream ss;
		ss << "# rank " << rank << " attempting to use GPU " << device << '\n';
		std::cout << ss.str();
	}
	if(m_gpuYield)
	{
		std::cout << "# setting CUDA to yield to other threads.\n";
		CUDA_SAFE_CALL( cudaSetDeviceFlags ( cudaDeviceScheduleYield ) );
	}
	return device;
}
