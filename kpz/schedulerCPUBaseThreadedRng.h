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

#ifndef KPZ_SCHEDULER_CPU_BASE_THREADED_RNG_H
#define KPZ_SCHEDULER_CPU_BASE_THREADED_RNG_H

#include "schedulerCPUBaseThreaded.h"

namespace Kpz
{
	template<class Rng>
	class SchedulerCPUBaseThreadedRng : public SchedulerCPUBaseThreaded
	{
		void init() {
			Rng::initialize(m_rngs);
		}
		protected:

		std::vector<Rng> m_rngs;

		/*! \tparam Scheduler must be derived from SchedulerCPUBaseThreadedRng<Rng> . */
		template<class Scheduler>
		struct Func_MCSBase {
			unsigned int m_id;
			Scheduler* m_this;
			Rng* m_rng;

			inline void pqSyncRngUpdate(int x, int y);
			inline void disorder2SyncRngUpdate(int x, int y);

				Func_MCSBase(int id, Scheduler* pthis);
				~Func_MCSBase();
		};
		
		public:
		
			SchedulerCPUBaseThreadedRng(int lx, int ly, unsigned int n_threads)
				: SchedulerCPUBaseThreaded(lx,ly,n_threads)
				  , m_rngs(nThreads())
			{init();}

			SchedulerCPUBaseThreadedRng(const SystemSize& size, unsigned int n_threads)
				: SchedulerCPUBaseThreaded(size,n_threads)
				  , m_rngs(nThreads())
			{init();}

		virtual void randomize() {Rng::randomize(m_rngs, m_dsfmt);}

		virtual bool collectH5(splash::DataCollector* data, const std::string& prefix = "");
		virtual bool retrieveH5(splash::DataCollector* data, const std::string& prefix = "");
	};
}

#include "schedulerCPUBaseThreadedRng.inl"
#endif
