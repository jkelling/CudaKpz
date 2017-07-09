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

#ifndef KMC_THREAD_POOL_H
#define KMC_THREAD_POOL_H

#include <vector>
#include <pthread.h>

namespace Kmc
{

class ThreadPool
{
	public:
	class ThreadData_base
	{
		unsigned int m_id;
		void* m_payload;
		pthread_t m_thread;
		ThreadPool* m_parent;

		void setPayload(void* payload) {m_payload = payload;}

		public:

			ThreadData_base(ThreadPool* parent, unsigned int id)
				: m_id(id), m_payload(0), m_thread(0), m_parent(parent) {}

		unsigned int id() const {return m_id;}
		void* const payload() const {return m_payload;}
		unsigned int nThreads() const {return m_parent->nThreads();}

		friend class ThreadPool;
	};

	private:
	std::vector<ThreadData_base> m_threads;
	unsigned int m_nThreads;

	static ThreadPool* s_global;
	static ThreadPool& initGLobalInstance();

		ThreadPool();
	public:

	template<class Payload>
	class ThreadData : public ThreadData_base
	{
			ThreadData();
		public:

		Payload* payload() {return reinterpret_cast<Payload*>(m_payload);}

		static ThreadData<Payload>* cast(void* p) {return reinterpret_cast<ThreadData<Payload>*>(p);}
	};

		ThreadPool(unsigned int n);
		~ThreadPool();

	void join();
	void setNThreadsMax(unsigned int n);
	unsigned int nThreadsMax() const {return m_threads.size();}
	unsigned int nThreads() const {return m_nThreads;}

	void run(void* payload, void* (*fct)(void*)) {run(payload, fct, nThreadsMax());}
	void run(void* payload, void* (*fct)(void*), unsigned int nThreads);

	static ThreadPool& global() {
		if(ThreadPool::s_global)
			return *ThreadPool::s_global;
		else
			return ThreadPool::initGLobalInstance();
	}
};

class ThreadInLinearDistribution
{
	size_t m_min, m_num;

	public:

		ThreadInLinearDistribution(const ThreadPool::ThreadData_base* thread, size_t n) {
			set(thread,n);
		}

	void set(const ThreadPool::ThreadData_base* thread, size_t n);

	inline size_t min() const {return m_min;}
	inline size_t num() const {return m_num;}
	inline size_t sup() const {return m_min+m_num;}
	inline size_t max() const {return sup()-1;}
};

}
#endif
