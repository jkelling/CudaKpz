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

#include "schedulerCPUpThread.h"

Kpz::SchedulerCPUpThread::SchedulerCPUpThread(int lx, int ly, int threads)
	: SchedulerCPU(lx, ly), m_numThreads((threads>0)?threads:1)
{
	init();
}
Kpz::SchedulerCPUpThread::SchedulerCPUpThread(const SystemSize& size, int threads)
	: SchedulerCPU(size), m_numThreads((threads>0)?threads:1)
{
	init();
}
void Kpz::SchedulerCPUpThread::init()
{
	if(m_numThreads <= 0)
		return;
	m_threads = new pthread_t[m_numThreads];
}

Kpz::SchedulerCPUpThread::~SchedulerCPUpThread()
{
	if(m_numThreads > 0)
	{
		joinMCS();
		delete m_threads;
	}
}

void Kpz::SchedulerCPUpThread::joinMCS()
{
	for(int a = 0; a < m_numThreads; ++a)
		pthread_join(m_threads[a], NULL);
}
