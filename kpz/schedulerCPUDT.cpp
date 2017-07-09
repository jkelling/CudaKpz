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

#ifndef KPZ_SCHEDULER_CPU_DT_CPP
#define KPZ_SCHEDULER_CPU_DT_CPP

#include "schedulerCPUDT.h"
#include "schedulerCPUDT_nested.h"

#include <ddIterators.h>
#include <shuffle.h>

template<class Rng>
const int Kpz::SchedulerCPUDT<Rng>::L_MIN[2] = {2,2};

template<class Rng>
Kpz::SchedulerCPUDT<Rng>::SchedulerCPUDT(int lx, int ly, int threads)
	: SchedulerCPUpThreadRng<Rng>(lx, ly, threads)
{
	std::cout << "SchedulerCPUDT\n";
	layoutFromL1();
}

template<class Rng>
Kpz::SchedulerCPUDT<Rng>::SchedulerCPUDT(const Kpz::SystemSize& size, int threads)
	: SchedulerCPUpThreadRng<Rng>(size, threads)
{
	std::cout << "SchedulerCPUDT\n";
	layoutFromL1();
}

template<class Rng>
void Kpz::SchedulerCPUDT<Rng>::layoutFromL1()
{
	const int size[2] = {SchedulerCPU::m_size.lDimX(), SchedulerCPU::m_size.lDimY()};
	m_dd.layoutFromL1(size, L_MIN, 1);
}

template<class Rng>
void Kpz::SchedulerCPUDT<Rng>::setLayout(const int layout[2])
{
	const int size[2] = {SchedulerCPU::m_size.lDimX(), SchedulerCPU::m_size.lDimY()};
	m_dd.setLayout(layout, size);
}

template<class Rng>
void Kpz::SchedulerCPUDT<Rng>::mcs(int n)
{
	// prepare functors for threads
	Func_MCSBase* func = (Func_MCSBase*)malloc(sizeof(Func_MCSBase)*numThreads());
	for(int a = 0; a < numThreads(); ++a)
		new (&func[a]) Func_MCSBase(a, this);
	
	// loop over mcs
	for(int a = 0; a < n; ++a)
	{
		Kmc::Shuffle<DD::N_DOMAIN_SET> seq(&m_dsfmt);

		for(int d = 0; d < DD::N_DOMAIN_SET; ++d)
		{
			m_dd.setDomainSet(seq[d]);
			for(int a = 0; a < numThreads(); ++a)
			{
				func[a].resetIter();
				run(&SchedulerCPUpThread::func_mcs<Func_MCSBase>, &func[a], a);
			}
			joinMCS();
		}
		++m_deviceMcs;
	}
	free(func);
}

template<class Rng>
std::ostream& Kpz::SchedulerCPUDT<Rng>::printLayout(std::ostream& out)
{
	out << "DT layout\n"
		"\tblockDim " << m_dd.blockDim(0) << 'x' << m_dd.blockDim(1) << "\n"
		"\tblocks " << (1<<m_dd.lBlocks(0)) << 'x' << (1<<m_dd.lBlocks(1))
		<< " =" << m_dd.nBlocks() << "\n"
		"\tRunning " << numThreads() << " threads per scheduler.\n";
			
	return out;
}

#endif
