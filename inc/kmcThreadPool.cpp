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

#include "kmcThreadPool.h"

#include <cstdlib>
#include <iostream>

#include <unistd.h>

Kmc::ThreadPool* Kmc::ThreadPool::s_global = 0;

Kmc::ThreadPool::~ThreadPool()
{
	join();
}

void Kmc::ThreadPool::join()
{
	for(unsigned int a = 0; a < m_nThreads; ++a)
	{
		pthread_join(m_threads[a].m_thread, NULL);
		m_threads[a].m_thread = 0;
	}
}

void Kmc::ThreadPool::setNThreadsMax(unsigned int n)
{
	join();
	if(n > m_threads.size())
	{
		m_threads.reserve(n);
		for(unsigned int a = m_threads.size(); a < n; ++a)
			m_threads.emplace_back(this, a);
	}
	else
		m_threads.erase(m_threads.begin()+n, m_threads.end());
	m_nThreads = n;
}

void Kmc::ThreadPool::run(void* payload, void* (*fct)(void*), unsigned int nThreads )
{
	join();
	m_nThreads = nThreads;
	for(unsigned int a = 0; a < m_nThreads; ++a)
	{
		m_threads[a].setPayload(payload);
		pthread_create(&m_threads[a].m_thread, NULL, fct, &m_threads[a]);
	}
}

Kmc::ThreadPool::ThreadPool(unsigned int n)
	: m_threads(), m_nThreads(0)
{
	setNThreadsMax(n);
}

Kmc::ThreadPool& Kmc::ThreadPool::initGLobalInstance()
{
	const char* env = std::getenv("OMP_NUM_THREADS");
	int n = 0;
	if(env)
		n = std::atoi(env);
	if(n <= 0)
		n = sysconf( _SC_NPROCESSORS_ONLN );
	ThreadPool::s_global = new ThreadPool(n);
	return *ThreadPool::s_global;
}

void Kmc::ThreadInLinearDistribution::set(const ThreadPool::ThreadData_base* thread, size_t n)
{
	const size_t div = n/thread->nThreads();
	const size_t mod = n%thread->nThreads();

	m_min = div*thread->id();
	m_num = div;
	if(mod > 0)
	{
		if(thread->id() < mod)
		{
			m_num += 1;
			m_min += thread->id();
		}
		else
		{
			m_min += mod;
		}
	}
}
