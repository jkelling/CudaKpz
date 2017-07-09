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

#ifndef KPZ_SCHEDULER_CPU_DB_NESTED_H
#define KPZ_SCHEDULER_CPU_DB_NESTED_H

#include "schedulerCPUDB.h"

#include <ddIterators.h>
#include <benchmarking.h>

namespace Kpz
{
	template<class Rng>
	class SchedulerCPUDB<Rng>::Func_MCSBase
	{
		protected:
		SchedulerCPUDB* m_this;
		int m_id;
		Kmc::DDIterator<DD> m_iter;

		int m_activeX;
		int m_activeY;
		int m_activeSites;

		public:
		Func_MCSBase(int id, SchedulerCPUDB* pthis)
			: m_this(pthis), m_id(id), m_iter(pthis->m_dd, id)
		{
			m_activeX = m_this->m_dd.blockDim(0)-DEAD_BORDER_DIM;
			m_activeY = m_this->m_dd.blockDim(1)-DEAD_BORDER_DIM;
			m_activeSites = m_activeX*m_activeY;
		}
		Func_MCSBase(int id, const Func_MCSBase& other)
			: m_this(other.m_this), m_id(id), m_iter(m_this->m_dd, id)
			, m_activeX(other.m_activeX), m_activeY(other.m_activeY), m_activeSites(other.m_activeSites)
		{}

		inline void exec() {
			// timeAinit;
			/*! \todo Why are the block updates slower by ~10% since template Rng is used? */
			do {
				// timeAstart;
				for(long long a = 0; a < m_activeSites; ++a)
				{
					// perform dead border update
					const int x = ((int)(m_this->rng(m_id).close_open()*m_activeX)
							+ m_iter.blockOffset(0))&m_this->m_sizeCache.mDimX();
					const int y = ((int)(m_this->rng(m_id).close_open()*m_activeY)
							+ m_iter.blockOffset(1))&m_this->m_sizeCache.mDimY();

					m_this->pqOneUpdate(x,y);
				}
				// timeAstopS("blockUpdate " << m_iter.index() << " thread " << m_id);
			} while(m_iter.next(m_this->numThreads()));
		}

		void resetIter() {
			m_iter.setIndex(m_id);
		}
	};
}

#endif
