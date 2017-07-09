/***************************************************************************
*   Copyright 2015 Jeffrey Kelling <j.kelling@hzdr.de>
*                  Helmholtz-Zentrum Dresden-Rossendorf
*                  Institute of Ion Beam Physics and Materials Research
*
*   This program is free software; you can redistribute it and/or
*   modify it under the terms of the GNU General Public
*   License as published by the Free Software Foundation; either
*   version 2 of the License, or (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program; if not, write to the Free Software
*   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
***************************************************************************/

#include <cstring>

template<class Rng>
	template<class Scheduler>
inline void Kpz::SchedulerCPUBaseThreadedRng<Rng>::Func_MCSBase<Scheduler>::pqSyncRngUpdate(int x, int y)
{
	LatticePoint c = m_this->m_size.latticePoint(x, y);
	Slope r = m_this->m_size.slopeX((x+1)&m_this->m_sizeCache.mDimX(), y);
	Slope l = m_this->m_size.slopeY(x, (y+1)&m_this->m_sizeCache.mDimY());

	const int config = c(m_this->m_stash) | (r(m_this->m_stash)<<2) | (l(m_this->m_stash)<<3);

	const float rnd = m_rng->randomFloat();
	if(config == 0b1100)
	{ // p
		if(rnd >= m_this->m_disorderP[0])
			return;
	}
	else if(config == 0b0011)
	{ // q
		if(rnd >= m_this->m_disorderQ[0])
			return;
	}
	else return;
	c.flip(m_this->m_stash);
	r.flip(m_this->m_stash);
	l.flip(m_this->m_stash);
}

template<class Rng>
	template<class Scheduler>
inline void Kpz::SchedulerCPUBaseThreadedRng<Rng>::Func_MCSBase<Scheduler>::disorder2SyncRngUpdate(int x, int y)
{
	LatticePoint c = m_this->m_size.latticePoint(x, y);
	Slope r = m_this->m_size.slopeX((x+1)&m_this->m_sizeCache.mDimX(), y);
	Slope l = m_this->m_size.slopeY(x, (y+1)&m_this->m_sizeCache.mDimY());

	const int config = c(m_this->m_stash) | (r(m_this->m_stash)<<2) | (l(m_this->m_stash)<<3);

	const float rnd = m_rng->randomFloat();
	if(config == 0b1100)
	{ // p
		if(rnd >= m_this->m_disorderP[c(m_this->m_disorder)])
			return;
	}
	else if(config == 0b0011)
	{ // q
		if(rnd >= m_this->m_disorderQ[c(m_this->m_disorder)])
			return;
	}
	else return;
	c.flip(m_this->m_stash);
	r.flip(m_this->m_stash);
	l.flip(m_this->m_stash);
}

template<class Rng>
	template<class Scheduler>
Kpz::SchedulerCPUBaseThreadedRng<Rng>::Func_MCSBase<Scheduler>::Func_MCSBase(int id, Scheduler* pthis)
	: m_id(id), m_this(pthis), m_rng(0)
{
}

template<class Rng>
	template<class Scheduler>
Kpz::SchedulerCPUBaseThreadedRng<Rng>::Func_MCSBase<Scheduler>::~Func_MCSBase()
{
	if(m_rng)
	{
		memcpy(&m_this->m_rngs[m_id], m_rng, sizeof(Rng));
		delete m_rng;
	}
}
