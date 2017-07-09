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

// #include "schedulerCPUBaseThreadedRng.h"

#include <KMCsplash.h>

template<class Rng>
bool Kpz::SchedulerCPUBaseThreadedRng<Rng>::collectH5(splash::DataCollector* data, const std::string& prefix)
{
	if(!SchedulerCPUBaseThreaded::collectH5(data, prefix))
		return false;
	Rng::writeH5(data, m_mcs, m_rngs, prefix);
	return true;
}

template<class Rng>
bool Kpz::SchedulerCPUBaseThreadedRng<Rng>::retrieveH5(splash::DataCollector* data, const std::string& prefix)
{
	if(!SchedulerCPUBaseThreaded::retrieveH5(data, prefix))
		return false;
	if(!Rng::readH5(data, m_mcs, m_rngs, prefix))
	{
		randomize();
		std::cout << "[SchedulerCPUBaseThreadedRng][retrieveH5][WW] Failed to restore GPU RNG state from file, randomizing. Check if dSFMT was restored.\n";
	}
	return true;
}
