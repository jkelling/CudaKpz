/***************************************************************************
*   Copyright 2013 Jeffrey Kelling <j.kelling@hzdr.de>
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

#include "schedulerServiceStash.h"

#include <cstring>

Kpz::SchedulerServiceStash::SchedulerServiceStash(int lx, int ly)
	: SchedulerService(lx, ly), m_stash(new unsigned int[m_size.sizeW()])
{
}
Kpz::SchedulerServiceStash::SchedulerServiceStash(const SystemSize& size)
	: SchedulerService(size), m_stash(new unsigned int[m_size.sizeW()])
{
}

void Kpz::SchedulerServiceStash::pushSystem()
{
	join();
	memcpy(m_stash, m_system, m_size.sizeW()<<2);
	m_deviceMcs = m_mcs;
}

void Kpz::SchedulerServiceStash::popSystem()
{
	join();
	memcpy(m_system, m_stash, m_size.sizeW()<<2);
	m_mcs = m_deviceMcs;
}

bool Kpz::SchedulerServiceStash::changedSizeEvent()
{
	if(!SchedulerService::changedSizeEvent())
		return false;
	delete m_stash;
	m_stash = new unsigned int[m_size.sizeW()];
	return true;
}
