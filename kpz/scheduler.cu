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

#include "scheduler.h"

#include "kpzConst.h"
#include "kpz.h"
#include "systemSize.h"
#include "LocalLayout.h"
#include "LocalLayoutIO.h"

#include <kmcExceptCUDA.h>

#include <utility>

#include <cuda.h>

template <class GpuRng, class LocalLayout>
void Kpz::Scheduler<GpuRng, LocalLayout>::doMcs(int n, Kpz::Scheduler<GpuRng, LocalLayout>::kernelCallerFct mcs)
{
	for(; n > 0; --n)
	{
		for(int a = 0; a < (1<<MyLocalLayout::L_MCS_DIV); ++a)
		{
			const int xw = (int)(dsfmt_genrand_close_open(&m_dsfmt)*(1<<MyLocalLayout::L_BLOCK_DIM_X));
			const int yw = (int)(dsfmt_genrand_close_open(&m_dsfmt)*(1<<MyLocalLayout::L_BLOCK_DIM_Y));
			(this->*mcs)(xw,yw);
		}
		++m_deviceMcs;
	}
	cudaThreadSynchronize();
	cudaError_t error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		throw KmcExcept::CUDAError(error, __FILE__, __LINE__);
	}
}

template <class GpuRng, class LocalLayout>
double Kpz::Scheduler<GpuRng, LocalLayout>::mcsScale() const
{
	const double sitesPerBlock = MyLocalLayout::BLOCK_DATA<<4;
	const double activeSitesPerBlock = (MyLocalLayout::BLOCK_BORDER_X)*
		(MyLocalLayout::BLOCK_BORDER_Y);
	const double scale = (activeSitesPerBlock)/sitesPerBlock;
	return scale;
}

template <class GpuRng, class LocalLayout>
void Kpz::Scheduler<GpuRng, LocalLayout>::init()
{
	std::cout << "Scheduler (GPU, DB4 ur)\n";
	MyLocalLayout::print(std::cout);
	m_size.setGPUfromLocalLayout<MyLocalLayout>();
	initStaticDeviceConstDyn(m_size());
	std::cout << "\nSet system size to " << m_size.lDimX() << ',' << m_size.lDimY() << '\n' << m_size;
}

template <class GpuRng, class LocalLayout>
Kpz::Scheduler<GpuRng, LocalLayout>::Scheduler(int lx, int ly, int device)
	: SchedulerBase(lx, ly, device), m_gpuRng(getRandomNumberGeneratorCount(MyLocalLayout::L_THREADS), &m_dsfmt)
	  , m_randomNumberBlockMultiplier(calcRandomNumberBlockMultiplier())
{
	init();
}

template <class GpuRng, class LocalLayout>
Kpz::Scheduler<GpuRng, LocalLayout>::Scheduler(const SystemSize& size, int device)
	: SchedulerBase(size, device), m_gpuRng(getRandomNumberGeneratorCount(MyLocalLayout::L_THREADS), &m_dsfmt)
	  , m_randomNumberBlockMultiplier(calcRandomNumberBlockMultiplier())
{
	init();
}

template <class GpuRng, class LocalLayout>
Kpz::Scheduler<GpuRng, LocalLayout>::Scheduler(SchedulerService &&other, int device)
	: SchedulerBase(std::move(other), device), m_gpuRng(getRandomNumberGeneratorCount(MyLocalLayout::L_THREADS), &m_dsfmt)
  , m_randomNumberBlockMultiplier(calcRandomNumberBlockMultiplier())
{
	init();
}

template <class GpuRng, class LocalLayout>
bool Kpz::Scheduler<GpuRng, LocalLayout>::changedSizeEvent()
{
	if(!SchedulerBase::changedSizeEvent())
		return false;
	m_randomNumberBlockMultiplier = calcRandomNumberBlockMultiplier();
	initStaticDeviceConstDyn(m_size());
	return true;
}

#include <KMCsplash.h>

template <class GpuRng, class LocalLayout>
bool Kpz::Scheduler<GpuRng, LocalLayout>::collectH5(splash::DataCollector* data, const std::string& prefix)
{
	if(!SchedulerBase::collectH5(data, prefix))
		return false;
	m_gpuRng.writeH5(data, m_mcs, prefix);
	return true;
}

template <class GpuRng, class LocalLayout>
bool Kpz::Scheduler<GpuRng, LocalLayout>::retrieveH5(splash::DataCollector* data, const std::string& prefix)
{
	if(!SchedulerBase::retrieveH5(data, prefix))
		return false;
	m_gpuRng.readH5(data, m_mcs, prefix);
	return true;
}

#include "schedulerInstances.h"
