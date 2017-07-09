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

#include "timer.h"

#include <cmath>

Kmc::Timer* Kmc::TimerSingleton::m_instance = 0;

void Kmc::Timer::set(float th)
{
	time_t now = time(0);
	tm* t = localtime(&now);
	const int h = th;
	th = (th-h)*60;
	const int min = th;
	th = (th-min)*60;
	const int sec = th;

	t->tm_hour += h;
	t->tm_min += min;
	t->tm_sec += sec;

	m_end = mktime(t);
}
