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

#pragma once

#include "schedulerBase.h"

namespace Kpz
{
	template<class GpuRng>
	class SchedulerSCABit : public SchedulerBase
	{
		protected:
		GpuRng* m_gpuRng;
		int m_updatesPerGen, m_lThreads;

		inline void setPQDefault() {
			setPQ(.5,.0);
		}

		typedef void (SchedulerSCABit::*kernelCallerFct) (int parity);
		void caller_mcsSyncRng(int parity);
		void caller_mcsPQSyncRng(int parity);
		template<class KpzCore> void run(int parity, KpzCore core);
		// void caller_mcsDisorder2SyncRng(int parity);

		void doMcs(int n, kernelCallerFct mcs);

		void initStaticDeviceConstDyn(const Kpz::SystemSize& size);
		void init();
		virtual bool changedSizeEvent();

		public:
			SchedulerSCABit(int lx, int ly, int device = -1);
			SchedulerSCABit(const SystemSize& size, int device = -1);
			SchedulerSCABit(SchedulerService &&other, int device = -1);

			virtual ~SchedulerSCABit();

		void randomize() {m_gpuRng->randomize();}
		virtual void initHomogeneous(Encoding encoding = ENC_CBBIT);

		void mcs(int n) {
			SchedulerSCABit::mcsSyncRng(n);
		}
		void mcsSyncRng(int n) {
			setPQDefault();
			doMcs(n, &SchedulerSCABit::caller_mcsSyncRng);
		}

		void mcsPQSyncRng(int n) {
			doMcs(n, &SchedulerSCABit::caller_mcsPQSyncRng);
		}
		void mcsDisorder2SyncRng(int n) {
			// doMcs(n, caller_mcsDisorder2SyncRng);
			std::cerr << "[BUG] mcsDisorder2SyncRng() not implemented.\n";
		}

		virtual double mcsScale() const;

		virtual bool collectH5(splash::DataCollector* data, const std::string& prefix = "");
		virtual bool retrieveH5(splash::DataCollector* data, const std::string& prefix = "");

		virtual void correlate();
	};
}
