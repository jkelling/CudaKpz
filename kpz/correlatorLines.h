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

#ifndef KPZ_CORRELATOR_LINES_H
#define KPZ_CORRELATOR_LINES_H

#include "correlator.h"
#include "schedulerService.h"

namespace Kpz
{
	class CorrelatorLines : public Correlator
	{
		double* m_hL;

			CorrelatorLines(const CorrelatorLines&); //no copy

		public:

			CorrelatorLines(SchedulerService* scheduler);
			virtual ~CorrelatorLines() {delete m_hL;}

		// void set(int mcs, Kpz::Roughness* r);
		// void correlate(Kpz::Roughness* r);
		virtual Kpz::Roughness set(int mcs);
		virtual Kpz::Roughness correlate();
	
	};
}

#endif
