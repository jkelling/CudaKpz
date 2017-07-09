/***************************************************************************
*   Copyright 2011 - 2014 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef KPZ_SCHEDULER_BASE_H
#define KPZ_SCHEDULER_BASE_H

#include "schedulerService.h"
#include "schedulerGPU.h"
#include "systemSize.h"

#include <utility>

namespace Kpz
{
	class SchedulerBase : public SchedulerService, public SchedulerGPU
	{
		protected:

		unsigned int *d_System;

		void initCUDA();

		static void* callMcs(void* args);

			SchedulerBase(int lx, int ly, int device, int smemWperBlock = 0, Encoding enc = ENC_LOCAL)
				: SchedulerService(lx, ly, enc), SchedulerGPU(device, smemWperBlock)
			{initCUDA(); m_callMcs = &callMcs;}
			SchedulerBase(const SystemSize& size, int device, int smemWperBlock = 0, Encoding enc = ENC_LOCAL)
				: SchedulerService(size, enc), SchedulerGPU(device, smemWperBlock)
			{initCUDA(); m_callMcs = &callMcs;}
			SchedulerBase(SchedulerService &&other, int device, int smemWperBlock = 0)
				: SchedulerService(std::move(other)), SchedulerGPU(device, smemWperBlock)
			{initCUDA(); m_callMcs = &callMcs;}
		public:
		
			~SchedulerBase();

		void pushSystem();
		void popSystem();

		virtual void mcsNop(int n) {m_deviceMcs += n;}

		virtual void setMcs(unsigned int mcs);

		/*!
		 * Currently we only really need the synched variant, not implmenting
		 * this one separately.
		 */
		virtual void mcsPQ(int n) {
			mcsPQSyncRng(n); 
		}
	};
}
#endif
