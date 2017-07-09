/***************************************************************************
*   Copyright 2011 - 2016 Jeffrey Kelling <j.kelling@hzdr.de>
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

#include "readWriteSemaphore.h"

ReadWriteSemaphore::ReadWriteSemaphore()
	: m_state(S_IDLE), m_locks(0)
{
	pthread_mutex_init(&m_mutex, NULL);
	pthread_cond_init(&m_idle, NULL);
	pthread_cond_init(&m_write, NULL);
}

ReadWriteSemaphore::~ReadWriteSemaphore()
{
	pthread_mutex_destroy(&m_mutex);
	pthread_cond_destroy(&m_idle);
	pthread_cond_destroy(&m_write);
}

void ReadWriteSemaphore::beginRead()
{
	pthread_mutex_lock(&m_mutex);
	if(m_state == S_WRITE)
		pthread_cond_wait(&m_idle, &m_mutex);
	++m_locks;
	m_state = S_READ;
	pthread_mutex_unlock(&m_mutex);
}

void ReadWriteSemaphore::beginWrite()
{
	pthread_mutex_lock(&m_mutex);
	switch(m_state)
	{
		case S_IDLE:
			m_state = S_WRITE;
			break;
		case S_READ:
			m_state = S_WRITE;
			pthread_cond_wait(&m_write, &m_mutex);
			m_locks = 1;
			break;
		case S_WRITE:
			++m_locks;
			pthread_cond_wait(&m_write, &m_mutex);
			break;
	}
	pthread_mutex_unlock(&m_mutex);
}

void ReadWriteSemaphore::endRead()
{
	pthread_mutex_lock(&m_mutex);
	--m_locks;
	if(m_locks == 0)
	{
		if(m_state == S_WRITE)
			pthread_cond_signal(&m_write);
		else
		{
			m_state = S_IDLE;
			pthread_cond_broadcast(&m_idle);
		}
	}
	pthread_mutex_unlock(&m_mutex);
}

void ReadWriteSemaphore::endWrite()
{
	pthread_mutex_lock(&m_mutex);
	--m_locks;
	if(m_locks == 0)
	{
		m_state = S_IDLE;
		pthread_cond_broadcast(&m_idle);
	}
	else
		pthread_cond_signal(&m_write);
	pthread_mutex_unlock(&m_mutex);
}

void ReadWriteSemaphore::wait()
{
	if(m_state != S_IDLE)
	{
		pthread_cond_wait(&m_idle, &m_mutex);
		pthread_mutex_unlock(&m_mutex);
	}
}
