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

#ifndef KPZ_SCHEDULER_CPU_DT_NESTED_H
#define KPZ_SCHEDULER_CPU_DT_NESTED_H

#include "schedulerCPUDT.h"

#include <ddIterators.h>
#include <benchmarking.h>

namespace Kpz
{
	template<class Rng>
	class SchedulerCPUDT<Rng>::Func_MCSBase
	{
		protected:
		SchedulerCPUDT<Rng>* m_this;
		int m_id;
		Kmc::DDIterator<DD> m_iter;

		public:
		Func_MCSBase(int id, SchedulerCPUDT<Rng>* pthis)
			: m_this(pthis), m_id(id), m_iter(pthis->m_dd, id)
		{
		}

		inline void exec() {
			// timeAinit;
			/*! \todo Why are the block updates slower by ~10% since template Rng is used? */
			const int blockSites = m_this->m_dd.blockSites();
			const int blockDimX = m_this->m_dd.blockDim(0);
			const int blockDimY = m_this->m_dd.blockDim(1);
			do {
				// timeAstart;
				for(long long a = 0; a < blockSites; ++a)
				{
					// perform dead border update
					const int x = ((int)(m_this->rng(m_id).close_open()*blockDimX)
							+ m_iter.blockOffset(0));
					const int y = ((int)(m_this->rng(m_id).close_open()*blockDimY)
							+ m_iter.blockOffset(1));

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
