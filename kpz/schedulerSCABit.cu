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
#include "kpz.h"
#include "systemSize.h"
#include "correlator.h"

#include <benchmarking.h>
#include <timer.h>
#include <kmcExceptCUDA.h>

#include <utility>
#include <sstream>
#include <iomanip>

#include <cuda.h>

template<class GpuRng>
double Kpz::SchedulerSCABit<GpuRng>::mcsScale() const
{
	return m_disorderP[0]+m_disorderQ[0];
}

template <class GpuRng>
void Kpz::SchedulerSCABit<GpuRng>::doMcs(int n, Kpz::SchedulerSCABit<GpuRng>::kernelCallerFct mcs)
{
	cudaEvent_t event;
	cudaEventCreate(&event);

	if(n <= 0)
		return;
	(this->*mcs)(0);
	(this->*mcs)(1);
	++m_deviceMcs;
	for(--n; n > 0; --n)
	{
		cudaEventRecord(event);
		(this->*mcs)(0);
		(this->*mcs)(1);
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
void Kpz::SchedulerSCABit<GpuRng>::init()
{
	std::cout << "SchedulerSCABit (GPU)\n";
	initStaticDeviceConstDyn(m_size);
	setPQDefault();
}

template <class GpuRng>
Kpz::SchedulerSCABit<GpuRng>::SchedulerSCABit(int lx, int ly, int device)
	: SchedulerBase(lx, ly, device, 0, ENC_CBBIT), m_gpuRng(0)
{
	init();
}

template <class GpuRng>
Kpz::SchedulerSCABit<GpuRng>::SchedulerSCABit(const SystemSize& size, int device)
	: SchedulerBase(size, device, 0, ENC_CBBIT), m_gpuRng(0)
{
	init();
}

template <class GpuRng>
Kpz::SchedulerSCABit<GpuRng>::SchedulerSCABit(SchedulerService &&other, int device)
	: SchedulerBase(std::move(other), device), m_gpuRng(0)
{
	init();
}

template <class GpuRng>
Kpz::SchedulerSCABit<GpuRng>::~SchedulerSCABit()
{
	delete m_gpuRng;
}

template <class GpuRng>
bool Kpz::SchedulerSCABit<GpuRng>::changedSizeEvent()
{
	if(!SchedulerBase::changedSizeEvent())
		return false;
	initStaticDeviceConstDyn(m_size);
	return true;
}

template <class GpuRng>
void Kpz::SchedulerSCABit<GpuRng>::correlate()
{
	if((!m_fAnalysisOnGPU) || m_mcs != m_deviceMcs)
		return SchedulerService::correlate();

	timeAinit;
	timeAstart;

	std::ostringstream text;
	bool corr = false;
	using ::operator<<;

	for(typename std::list<Correlator*>::iterator a = m_correlateTags.begin(); a != m_correlateTags.end();)
	{
		if(m_mcs > (*a)->stop())
			a = m_correlateTags.erase(a);
		else
		{
			joinMcs();
			auto m_roughness = (*a)->correlateCuda(d_System, m_blocks);
			if(!m_silentRoughness)
				text << std::setprecision(m_outPrecision) << m_mcs << '\t' << m_roughness << '\t' << (**a) << '\n';
			++a;
			corr = true;
		}
	}

	if(!corr)
	{
		return roughnessAsync();
	}
	if(!m_silentRoughness)
	{
		cout() << text.str();
		cout().flush();
	}

	timeAstop("SchedulerSCABit::correlate(correlateCuda)");
}

#include <KMCsplash.h>

template <class GpuRng>
bool Kpz::SchedulerSCABit<GpuRng>::collectH5(splash::DataCollector* data, const std::string& prefix)
{
	if(!SchedulerBase::collectH5(data, prefix))
		return false;
	m_gpuRng->writeH5(data, m_mcs, prefix);
	return true;
}

template <class GpuRng>
bool Kpz::SchedulerSCABit<GpuRng>::retrieveH5(splash::DataCollector* data, const std::string& prefix)
{
	if(!SchedulerBase::retrieveH5(data, prefix))
		return false;
	if(!m_gpuRng->readH5(data, m_mcs, prefix))
	{
		m_gpuRng->randomize();
		std::cout << "[SchedulerSCABit][retrieveH5][WW] Failed to restore GPU RNG state from file, randomizing. Check if dSFMT was restored.\n";
	}
	return true;
}

template<class Rng>
void Kpz::SchedulerSCABit<Rng>::initHomogeneous(Encoding encoding)
{
	SchedulerService::initHomogeneous(ENC_CBBIT);
}

#include "schedulerSCABitInstances.h"
