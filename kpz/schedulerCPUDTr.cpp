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

#ifndef KPZ_SCHEDULER_CPU_DTR_CPP
#define KPZ_SCHEDULER_CPU_DTR_CPP

#include "schedulerCPUDTr.h"
#include "schedulerCPUDT_nested.h"

#include <shuffle.h>
#include <ddIterators.h>

template<class Rng>
const int Kpz::SchedulerCPUDTr<Rng>::L_MIN[2] = {2,2};

template<class Rng>
Kpz::SchedulerCPUDTr<Rng>::SchedulerCPUDTr(int lx, int ly, int threads)
	: SchedulerCPUBaseThreadedRng<Rng>(lx, ly, threads)
{
	std::cout << "SchedulerCPUDTr\n";
	layoutFromL1();
}

template<class Rng>
Kpz::SchedulerCPUDTr<Rng>::SchedulerCPUDTr(const Kpz::SystemSize& size, int threads)
	: SchedulerCPUBaseThreadedRng<Rng>(size, threads)
{
	std::cout << "SchedulerCPUDTr\n";
	layoutFromL1();
}

template<class Rng>
void Kpz::SchedulerCPUDTr<Rng>::layoutFromL1()
{
	const int size[2] = {SchedulerService::m_size.lDimX(), SchedulerService::m_size.lDimY()};
	m_dd.layoutFromL1(size, L_MIN, 1);
}

template<class Rng>
void Kpz::SchedulerCPUDTr<Rng>::setLayout(const int layout[2])
{
	const int size[2] = {SchedulerService::m_size.lDimX(), SchedulerService::m_size.lDimY()};
	m_dd.setLayout(layout, size);
}

template<class Rng>
	template<class UpdateFkt>
void Kpz::SchedulerCPUDTr<Rng>::doMcs(int n, UpdateFkt updateFkt)
{
	// prepare functors for threads
	Func_MCS<UpdateFkt>* func = (Func_MCS<UpdateFkt>*)malloc(sizeof(Func_MCS<UpdateFkt>)*nThreads());
	for(int a = 0; a < nThreads(); ++a)
		new (&func[a]) Func_MCS<UpdateFkt>(a, this, updateFkt);
	
	// loop over mcs
	for(int a = 0; a < n; ++a)
	{
		Kmc::Shuffle<DD::N_DOMAIN_SET> seq(&m_dsfmt);
		{
			const unsigned int shift[2] =
				{(unsigned int)(dsfmt_genrand_close_open(&m_dsfmt)*m_dd.blockDim(0))
					, (unsigned int)(dsfmt_genrand_close_open(&m_dsfmt)*m_dd.blockDim(1))};
			m_dd.setShift(shift);
		}

		for(int d = 0; d < DD::N_DOMAIN_SET; ++d)
		{
			m_dd.setDomainSet(seq[d]);
			for(int a = 0; a < nThreads(); ++a)
			{
				func[nextThread()].resetIter();
				run(&SchedulerCPUBaseThreaded::func_mcs<Func_MCS<UpdateFkt> >, &func[nextThread()]);

			}
			joinAllWorkers();
		}
		++m_deviceMcs;
	}
	free(func);
}

template<class Rng>
void Kpz::SchedulerCPUDTr<Rng>::mcsSyncRng(int n)
{
	joinRoughness();
	doMcs(n, [](MyFunc_MCSBase* f, int x, int y){f->m_this->pqOneUpdate(x,y);});
}

template<class Rng>
void Kpz::SchedulerCPUDTr<Rng>::mcsPQSyncRng(int n)
{
	joinRoughness();
	doMcs(n, [](MyFunc_MCSBase* f, int x, int y){f->pqSyncRngUpdate(x,y);});
}

template<class Rng>
void Kpz::SchedulerCPUDTr<Rng>::mcsDisorder2SyncRng(int n)
{
	joinRoughness();
	doMcs(n, [](MyFunc_MCSBase* f, int x, int y){f->disorder2SyncRngUpdate(x,y);});
}

template<class Rng>
std::ostream& Kpz::SchedulerCPUDTr<Rng>::printLayout(std::ostream& out)
{
	out << "DT layout\n"
		"\tblockDim " << m_dd.blockDim(0) << 'x' << m_dd.blockDim(1) << "\n"
		"\tblocks " << (1<<m_dd.lBlocks(0)) << 'x' << (1<<m_dd.lBlocks(1))
		<< " =" << m_dd.nBlocks() << "\n"
		"\tRunning " << nThreads() << " threads per scheduler.\n";
			
	return out;
}

template<class Rng>
	template<class UpdateFkt>
struct Kpz::SchedulerCPUDTr<Rng>::Func_MCS : public MyFunc_MCSBase
{
	using MyFunc_MCSBase::m_this;
	using MyFunc_MCSBase::m_rng;
	using MyFunc_MCSBase::m_id;

	UpdateFkt m_updateFkt;
	Kmc::DDIterator<DD> m_iter;

	Func_MCS(unsigned int id, SchedulerCPUDTr<Rng>* pthis, UpdateFkt updateFkt)
		: MyFunc_MCSBase(id, pthis)
		  , m_updateFkt(updateFkt), m_iter(m_this->m_dd)
	{}

	inline void exec() {
		Rng rng;
		memcpy(&rng, &m_this->m_rngs[m_id], sizeof(Rng));
		m_rng = &rng;
		// timeAinit;
		const int blockSites = m_this->m_dd.blockSites();
		const int blockDimX = m_this->m_dd.blockDim(0);
		const int blockDimY = m_this->m_dd.blockDim(1);
		// timeAstart;

		do {
			// timeAstart;
			for(long long a = 0; a < blockSites; ++a)
			{
				// perform dead border update
				const int x = ((int)(m_rng->randomFloat()*blockDimX)
						+ m_iter.blockOffset(0))&m_this->m_sizeCache.mDimX();
				const int y = ((int)(m_rng->randomFloat()*blockDimY)
						+ m_iter.blockOffset(1))&m_this->m_sizeCache.mDimX();

				m_updateFkt(this,x,y);
			}
			// timeAstopS("blockUpdate " << m_iter.index() << " thread " << m_id);
		} while(m_iter.next(m_this->nThreads()));

		memcpy(&m_this->m_rngs[m_id], m_rng, sizeof(Rng));
	}

	void resetIter() {
		m_iter.setIndex(m_id);
	}
};
#endif
