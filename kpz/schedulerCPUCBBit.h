/***************************************************************************
*   Copyright 2015 - 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#pragma once

#include "schedulerCPUCB.h"

namespace Kpz
{
	/*! CPU scheduler with lattce checkerboard DD (SCA) bit parallel updates.
	 *
	 * For now limited to p=.5
	 *
	 * No stash. Will join analysis before beginning MCS.
	 *
	 */
	template<class Rng>
	class SchedulerCPUCBBit : public SchedulerCPUCB<Rng>
	{
		protected:

		using SchedulerService::m_deviceMcs;
		using SchedulerService::m_dsfmt;
		using SchedulerService::m_size;
		using SchedulerService::m_disorderP;
		using SchedulerService::m_disorderQ;
		using SchedulerService::joinRoughness;
		using SchedulerService::mcsNop;
		using SchedulerService::initHomogeneous;
		using typename SchedulerService::Encoding;
		using SchedulerService::ENC_CBBIT;
		using SchedulerCPUBaseThreaded::run;
		using SchedulerCPUBaseThreaded::joinWorker;
		using SchedulerCPUBaseThreaded::joinAllWorkers;
		using SchedulerCPUBaseThreaded::nextThread;
		using SchedulerCPUCB<Rng>::m_signum;

		struct Func_MCSLocalBase;
		template<class UpdateFkt>
		struct Func_MCS;
		typedef typename SchedulerCPUBaseThreadedRng<Rng>::template Func_MCSBase<Kpz::SchedulerCPUCBBit<Rng> > MyFunc_MCSBase;

		template<class UpdateFkt>
		void doMcs(int n, UpdateFkt updateFkt);

		public:
		
			SchedulerCPUCBBit(int lx, int ly, int threads = 1);
			SchedulerCPUCBBit(const SystemSize& size, int threads = 1);

		virtual void initHomogeneous(Encoding encoding = ENC_CBBIT);

		virtual void mcs(int n);
		virtual void mcsSyncRng(int n);
		virtual void mcsPQSyncRng(int n);
		/*! \todo requires correct disorder encoding, which is not implemented
		 * in SchedulerService yet.
		 */
		virtual void mcsDisorder2SyncRng(int n);

		using SchedulerCPUBaseThreaded::nThreads;
	};

}
