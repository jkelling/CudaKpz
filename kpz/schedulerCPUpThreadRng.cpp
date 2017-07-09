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


template<class Rng>
Kpz::SchedulerCPUpThreadRng<Rng>::SchedulerCPUpThreadRng(int lx, int ly, int threads)
	: SchedulerCPUpThread(lx, ly, threads)
{
	m_rng = Rng::create(numThreads(), &m_dsfmt);
}
template<class Rng>
Kpz::SchedulerCPUpThreadRng<Rng>::SchedulerCPUpThreadRng(const SystemSize& size, int threads)
	: SchedulerCPUpThread(size, threads)
{
	m_rng = Rng::create(numThreads(), &m_dsfmt);
}
template<class Rng>
Kpz::SchedulerCPUpThreadRng<Rng>::~SchedulerCPUpThreadRng()
{
	delete m_rng;
}
