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

#ifndef KMC_READ_WRITE_SEMAPHORE_H
#define KMC_READ_WRITE_SEMAPHORE_H

#include <pthread.h>

class ReadWriteSemaphore
{
	public:
	enum State { S_IDLE, S_READ, S_WRITE };

	private:
	int m_locks;
	State m_state;
	pthread_mutex_t m_mutex;
	pthread_cond_t m_idle, m_write;

	public:

		ReadWriteSemaphore();
		~ReadWriteSemaphore();

	void beginRead();
	void beginWrite();
	void endRead();
	void endWrite();
	void wait();

	inline State state() const {return m_state;}
};

#endif
