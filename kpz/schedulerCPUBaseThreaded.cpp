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

#include "schedulerCPUBaseThreaded.h"

void Kpz::SchedulerCPUBaseThreaded::run(void* fkt(void*), void* arg)
{
	joinWorker(m_nextThread);
	pthread_create(&m_threads[m_nextThread], NULL, fkt, arg);
	++m_nextThread;
	if(m_nextThread>=m_threads.size())
		m_nextThread = 0;
}

unsigned int Kpz::SchedulerCPUBaseThreaded::joinAllWorkers()
{
	for(unsigned int a = m_nextThread; a < m_threads.size(); ++a)
		joinWorker(a);
	for(unsigned int a = 0; a < m_nextThread; ++a)
		joinWorker(a);
	m_nextThread = 0;
}

Kpz::SchedulerCPUBaseThreaded::SchedulerCPUBaseThreaded(int lx, int ly, unsigned int nThreads)
	: SchedulerServiceStash(lx,ly), m_sizeCache(m_size)
	  , m_threads(0), m_nextThread(0)
{
	m_threads.resize((nThreads>0)?nThreads:1,0);
}

Kpz::SchedulerCPUBaseThreaded::SchedulerCPUBaseThreaded(const SystemSize& size, unsigned int nThreads)
	: SchedulerServiceStash(size), m_sizeCache(m_size)
	  , m_threads(0), m_nextThread(0)
{
	m_threads.resize((nThreads>0)?nThreads:1,0);
}
