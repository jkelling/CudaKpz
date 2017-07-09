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

#ifndef KPZ_SCHEDULER_CPU_BASE_THREADED_H
#define KPZ_SCHEDULER_CPU_BASE_THREADED_H

#include "schedulerServiceStash.h"
#include "systemSize.h"
#include "param.h"

#include <benchmarking.h>

#include <vector>

#include <pthread.h>

#define ASYNC_BENCH(a) \
	timeAinit; \
	timeAstart; \
	a; \
	timeAstopS("SchedulerCPUBaseThreaded::mcs_" << n); \


namespace Kpz
{
	class SchedulerCPUBaseThreaded : public SchedulerServiceStash
	{
		protected:

		SystemSizeCPU m_sizeCache;
		std::vector<pthread_t> m_threads;
		unsigned int m_nextThread;

		void joinWorker(unsigned int t) {
			if(m_threads[t] != 0)
			{
				pthread_join(m_threads[t], NULL);
				m_threads[t] = 0;
			}
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

		// functions to be used for scheduling treads by derived class:
		unsigned int nextThread() const {return m_nextThread;}
		void run(void* fkt(void*), void* arg);
		unsigned int joinAllWorkers();


		inline void pqOneUpdate(int x, int y);

		public:
		
			SchedulerCPUBaseThreaded(int lx, int ly, unsigned int nThreads);
			SchedulerCPUBaseThreaded(const SystemSize& size, unsigned int nThreads);

		unsigned int nThreads() const {return m_threads.size();}

		virtual void mcs(int n) {mcsSyncRng(n);}
		virtual void mcsPQ(int n) {mcsPQSyncRng(n);}
		virtual void mcsPQAsync(int n) {mcsPQSyncRng(n);}
		virtual void mcsDisorder2(int n) {mcsDisorder2SyncRng(n);}
		virtual void mcsDisorder2Async(int n) {ASYNC_BENCH(mcsDisorder2SyncRng(n));}
		virtual void mcsSyncRngAsync(int n) {ASYNC_BENCH(mcsSyncRng(n));}
		virtual void mcsPQSyncRngAsync(int n) {ASYNC_BENCH(mcsPQSyncRng(n));}
		virtual void mcsDisorder2SyncRngAsync(int n) {ASYNC_BENCH(mcsDisorder2SyncRng(n));}
	};
}
#include "schedulerCPUBaseThreaded.inl"
#endif
