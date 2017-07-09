/***************************************************************************
*   Copyright 2013 - 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef KPZ_SCHEDULER_SERVICE_STASH_H
#define KPZ_SCHEDULER_SERVICE_STASH_H

#include "schedulerService.h"
#include "systemSize.h"

namespace Kpz
{
	class SchedulerServiceStash : public SchedulerService
	{
		protected:

		unsigned int* m_stash;

		virtual bool changedSizeEvent();
		
		public:
		
			SchedulerServiceStash(int lx, int ly);
			SchedulerServiceStash(const SystemSize& size);

		unsigned int* stash() {return m_stash;}

		void pushSystem();
		void popSystem();
	};

}
#endif
