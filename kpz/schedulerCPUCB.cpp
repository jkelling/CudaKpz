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

#ifndef KPZ_SCHEDULER_CPU_CB_CPP
#define KPZ_SCHEDULER_CPU_CB_CPP

#include "schedulerCPUCB.h"

template<class Rng>
Kpz::SchedulerCPUCB<Rng>::SchedulerCPUCB(int lx, int ly, int threads)
	: SchedulerCPUBaseThreadedRng<Rng>(lx, ly, threads)
{
	std::cout << "SchedulerCPUCB\n";
}

template<class Rng>
Kpz::SchedulerCPUCB<Rng>::SchedulerCPUCB(const Kpz::SystemSize& size, int threads)
	: SchedulerCPUBaseThreadedRng<Rng>(size, threads)
{
	std::cout << "SchedulerCPUCB\n";
}

template<class Rng>
	template<class UpdateFkt>
void Kpz::SchedulerCPUCB<Rng>::doMcs(int n, UpdateFkt updateFkt)
{
	// prepare functors for threads
	Func_MCS<UpdateFkt>* func = (Func_MCS<UpdateFkt>*)malloc(sizeof(Func_MCS<UpdateFkt>)*nThreads());
	{
		const int xSitesThread = (this->m_size.dimXW()/nThreads())<<Kpz::L_BASIC_CELL_DIM_X;
		const int xRemainderThreads = (this->m_size.dimXW()%nThreads());
		int xSup = 0;
		for(int a = 0; a < xRemainderThreads; ++a)
		{
			const int xMin = xSup;
			xSup += xSitesThread + Kpz::BASIC_CELL_DIM_X;
			new (&func[a]) Func_MCS<UpdateFkt>(a, xMin, xSup, this, updateFkt, &func[(a+1)]);
		}
		for(int a = xRemainderThreads; a < nThreads(); ++a)
		{
			const int xMin = xSup;
			xSup += xSitesThread;
			new (&func[a]) Func_MCS<UpdateFkt>(a, xMin, xSup, this, updateFkt, &func[(a+1)%nThreads()]);
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
void Kpz::SchedulerCPUCB<Rng>::mcs(int n)
{
	SchedulerService::setPQ(.95, .0);
	mcsPQSyncRng(n);
}

template<class Rng>
void Kpz::SchedulerCPUCB<Rng>::mcsSyncRng(int n)
{
	SchedulerCPUCB<Rng>::mcsPQSyncRng(n);
}

template<class Rng>
void Kpz::SchedulerCPUCB<Rng>::mcsPQSyncRng(int n)
{
	joinRoughness();
	doMcs(n, [](MyFunc_MCSBase* f, int x, int y){f->pqSyncRngUpdate(x,y);});
	// doMcs(n, &MyFunc_MCSBase::pqSyncRngUpdate); // passing member function directly is somehow slower on it16
}

template<class Rng>
void Kpz::SchedulerCPUCB<Rng>::mcsDisorder2SyncRng(int n)
{
	joinRoughness();
	doMcs(n, [](MyFunc_MCSBase* f, int x, int y){f->disorder2SyncRngUpdate(x,y);});
	// doMcs(n, &MyFunc_MCSBase::disorder2SyncRngUpdate);
}

#include <benchmarking.h>
#include <algorithm>
#include <mutex>

template<class Rng>
	template<class UpdateFkt>
class Kpz::SchedulerCPUCB<Rng>::Func_MCS : public MyFunc_MCSBase
{
	protected:
	using MyFunc_MCSBase::m_this;
	using MyFunc_MCSBase::m_rng;
	using MyFunc_MCSBase::m_id;

	int m_xMin, m_xSup;
	UpdateFkt m_updateFkt;
	Func_MCS* m_right;
	std::mutex m_leftMutex;

	void updateLoop(int xMin, int xSup) {
		const register int signum = m_this->m_signum;
		for(int yW = 0; yW < m_this->m_size.dimY(); yW+=4)
			for(int xW = xMin; xW < xSup; xW+=4)
			{
				for(int y = yW; y < yW+4; ++y)
				{
					int x = xW + (y&1)^signum;
					m_updateFkt(this,x,y);
					// (this->*m_updateFkt)(x,y);
					x += 2;
					m_updateFkt(this,x,y);
					// (this->*m_updateFkt)(x,y);
				}
			}
	}

	public:
	Func_MCS(unsigned int id, int xMin, int xSup, SchedulerCPUCB<Rng>* pthis, UpdateFkt updateFkt, Func_MCS* right = 0)
		: MyFunc_MCSBase(id, pthis), m_updateFkt(updateFkt)
		, m_xMin(xMin), m_xSup(xSup), m_right(right)
	{
	}

	void setRight(Func_MCS* right) {m_right(right);}

	inline void exec() {
		// timeAinit;
		Rng rng;
		memcpy(&rng, &m_this->m_rngs[m_id], sizeof(Rng));
		m_rng = &rng;
		// timeAstart;

		m_leftMutex.lock();
		const int xCenter = m_xMin + (((m_xSup-m_xMin)/2)&(~3)); // round to multiple of 4
		updateLoop(m_xMin, xCenter);
		m_leftMutex.unlock();
		m_right->m_leftMutex.lock();
		updateLoop(xCenter, m_xSup);
		m_right->m_leftMutex.unlock();

		// timeAstopS("updates thread (xMin): " << m_xMin);
		memcpy(&m_this->m_rngs[m_id], m_rng, sizeof(Rng));
	}
};
#endif
