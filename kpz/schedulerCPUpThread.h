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

#ifndef KPZ_SCHEDULER_CPU_P_THREAD_H
#define KPZ_SCHEDULER_CPU_P_THREAD_H

#include "schedulerCPU.h"
#include "systemSize.h"

#include <pthread.h>

namespace Kpz
{

	class SchedulerCPUpThread : public SchedulerCPU
	{
		int m_numThreads;
		pthread_t* m_threads;
		void init();

		protected:
		
			SchedulerCPUpThread(int lx, int ly, int threads = 1);
			SchedulerCPUpThread(const SystemSize& size, int threads = 1);

		void run(void* fkt(void*), void* arg, int id) {
			pthread_create(&m_threads[id], NULL, fkt, arg);
		}
		/*! Generic funktor for run(...). Wraps given predicate which is
		 * expeced to provide an function exec():
		 */
		template<class Pred>
		static void* func_mcs(void* argOfTypePred)
		{
			Pred* arg = reinterpret_cast<Pred*>(argOfTypePred);
			arg->exec();
		}


		public:
		
			~SchedulerCPUpThread();

		int numThreads() const {return m_numThreads;}
		void joinMCS();
	};

	template<class Rng>
	class SchedulerCPUpThreadRng : public SchedulerCPUpThread
	{
		Rng* m_rng;

		protected:
		
		Rng& rng(int i) {return m_rng[i];}
		
		public:
			SchedulerCPUpThreadRng(int lx, int ly, int threads = 1);
			SchedulerCPUpThreadRng(const SystemSize& size, int threads = 1);
			~SchedulerCPUpThreadRng();

		//! \todo rng sync
	};

}

#include "schedulerCPUpThreadRng.cpp"

#endif
