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

#include <benchmarking.h>
#include <algorithm>
#include <mutex>

template<class Rng>
class Kpz::SchedulerCPUCBBit<Rng>::Func_MCSLocalBase : public MyFunc_MCSBase
{
	protected:
	using MyFunc_MCSBase::m_this;
	using MyFunc_MCSBase::m_rng;
	using MyFunc_MCSBase::m_id;

	int m_yMin, m_ySup;
	Func_MCSLocalBase* m_next;
	std::mutex m_topMutex;

	public:
	Func_MCSLocalBase(unsigned int id, int yMin, int ySup, SchedulerCPUCBBit<Rng>* pthis, Func_MCSLocalBase* next = 0)
		: MyFunc_MCSBase(id, pthis)
		, m_yMin(yMin), m_ySup(ySup), m_next(next)
	{
	}

	void setNext(Func_MCSLocalBase* next) {m_next = next;}
	std::mutex& topMutex() {return m_topMutex;}

	void pSyncRngUpdateVL5(int yMin, int ySup);
	void pSyncRngUpdateVL4(int yMin, int ySup);
	void pSyncRngUpdateVL3(int yMin, int ySup);
	void pqSyncRngUpdate(int yMin, int ySup);
	void disorder2SyncRngUpdate(int yMin, int ySup);
};

template<class Rng>
	template<class UpdateFkt>
class Kpz::SchedulerCPUCBBit<Rng>::Func_MCS : public Func_MCSLocalBase
{
	protected:
	using MyFunc_MCSBase::m_this;
	using MyFunc_MCSBase::m_rng;
	using Func_MCSLocalBase::m_id;
	using Func_MCSLocalBase::m_topMutex;
	using Func_MCSLocalBase::m_ySup;
	using Func_MCSLocalBase::m_yMin;
	using Func_MCSLocalBase::m_next;

	UpdateFkt m_updateFkt;

	public:
	Func_MCS(unsigned int id, int yMin, int ySup, SchedulerCPUCBBit<Rng>* pthis, UpdateFkt updateFkt, Func_MCS* next = 0)
		: Func_MCSLocalBase(id, yMin, ySup, pthis, next), m_updateFkt(updateFkt)
	{
	}

	void exec() {
		// timeAinit;
		Rng rng;
		memcpy(&rng, &m_this->m_rngs[m_id], sizeof(Rng));
		m_rng = &rng;
		// timeAstart;

		m_topMutex.lock();
		const int yCenter = m_yMin + (m_ySup-m_yMin)/2;
		m_updateFkt(this, m_yMin, yCenter);
		m_topMutex.unlock();
		m_next->topMutex().lock();
		m_updateFkt(this, yCenter, m_ySup);
		m_next->topMutex().unlock();

		// timeAstopS("updates thread (xMin): " << m_xMin);
		memcpy(&m_this->m_rngs[m_id], m_rng, sizeof(Rng));
	}
};

#include "schedulerCPUCBBit_nested.inl"
