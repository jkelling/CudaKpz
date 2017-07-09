/***************************************************************************
*   Copyright 2015 - 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#pragma once

#include "schedulerCPUCBBit.h"
#include "schedulerCPUCBBit_nested.h"

template<class Rng>
Kpz::SchedulerCPUCBBit<Rng>::SchedulerCPUCBBit(int lx, int ly, int threads)
	: SchedulerCPUCB<Rng>(lx, ly, threads)
{
	std::cout << "SchedulerCPUCBBit\n";
}

template<class Rng>
Kpz::SchedulerCPUCBBit<Rng>::SchedulerCPUCBBit(const Kpz::SystemSize& size, int threads)
	: SchedulerCPUCB<Rng>(size, threads)
{
	std::cout << "SchedulerCPUCBBit\n";
}

template<class Rng>
	template<class UpdateFkt>
void Kpz::SchedulerCPUCBBit<Rng>::doMcs(int n, UpdateFkt updateFkt)
{
	// prepare functors for threads
	Func_MCS<UpdateFkt>* func = (Func_MCS<UpdateFkt>*)malloc(sizeof(Func_MCS<UpdateFkt>)*nThreads());
	{
		const int ySitesThread = (this->m_size.dimY()/nThreads());
		const int yRemainderThreads = (this->m_size.dimY()%nThreads());
		int ySup = 0;
		for(int a = 0; a < yRemainderThreads; ++a)
		{
			const int yMin = ySup;
			ySup += ySitesThread + 1;
			new (&func[a]) Func_MCS<UpdateFkt>(a, yMin, ySup, this, updateFkt, &func[(a+1)]);
		}
		for(int a = yRemainderThreads; a < nThreads(); ++a)
		{
			const int yMin = ySup;
			ySup += ySitesThread;
			new (&func[a]) Func_MCS<UpdateFkt>(a, yMin, ySup, this, updateFkt, &func[(a+1)%nThreads()]);
		}
	}
	
	// loop over mcs
	for(int a = 0; a < n; ++a)
	{
		for(m_signum = 0; m_signum < 2; ++m_signum)
		{
			for(int a = 0; a < nThreads(); ++a)
			{
				run(&SchedulerCPUBaseThreaded::func_mcs<Func_MCS<UpdateFkt> >, &func[nextThread()]);
			}
			joinAllWorkers();
		}
		++m_deviceMcs;
	}
	free(func);
}

template<class Rng>
void Kpz::SchedulerCPUCBBit<Rng>::mcs(int n)
{
	mcsSyncRng(n);
}

template<class Rng>
void Kpz::SchedulerCPUCBBit<Rng>::mcsSyncRng(int n)
{
	joinRoughness();
	doMcs(n, [](Func_MCSLocalBase* f, int yMin, int ySup){f->pSyncRngUpdateVL3(yMin,ySup);});
}

template<class Rng>
void Kpz::SchedulerCPUCBBit<Rng>::mcsPQSyncRng(int n)
{
	joinRoughness();
	if(m_disorderQ[0] == 0 && m_disorderP[0] == .5)
		mcsSyncRng(n);
	else
		doMcs(n, [](Func_MCSLocalBase* f, int yMin, int ySup){f->pqSyncRngUpdate(yMin,ySup);});
}

template<class Rng>
void Kpz::SchedulerCPUCBBit<Rng>::mcsDisorder2SyncRng(int n)
{
	joinRoughness();
	doMcs(n, [](Func_MCSLocalBase* f, int yMin, int ySup){f->disorder2SyncRngUpdate(yMin,ySup);});
}

template<class Rng>
void Kpz::SchedulerCPUCBBit<Rng>::initHomogeneous(Encoding encoding)
{
	SchedulerService::initHomogeneous(ENC_CBBIT);
}
