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

#ifndef KPZ_SCHEDULER_SCA_H
#define KPZ_SCHEDULER_SCA_H

#include "schedulerBase.h"

namespace Kpz
{
	template<class GpuRng>
	class SchedulerSCA : public SchedulerBase
	{
		protected:
		GpuRng* m_gpuRng;
		int m_updatesPerGen;

		int m_lThreads, m_blocksX; // m_blocks; inherited from SchedulerGPU
		inline void setPQDefault() {
			setPQ(.95,0.);
		}

		typedef void (SchedulerSCA::*kernelCallerFct) (int parity);
		void caller_mcsPQSyncRng(int parity);
		void caller_mcsDisorder2SyncRng(int parity);

		void doMcs(int n, kernelCallerFct mcs);

		void initStaticDeviceConstDyn(const Kpz::SystemSize& size);
		void init();
		virtual bool changedSizeEvent();

		unsigned int* d_Disorder;
		void initDisorder2DeviceConst();

		public:
			SchedulerSCA(int lx, int ly, int device = -1);
			SchedulerSCA(const SystemSize& size, int device = -1);
			SchedulerSCA(SchedulerService &&other, int device = -1);

			~SchedulerSCA();

		void releaseDeviceDisorder();

		void randomize() {m_gpuRng->randomize();}

		void mcs(int n) {
			SchedulerSCA::mcsSyncRng(n);
		}
		void mcsSyncRng(int n) {
			setPQDefault();
			mcsPQSyncRng(n);
		}

		void mcsPQSyncRng(int n) {
			doMcs(n, &SchedulerSCA::caller_mcsPQSyncRng);
		}
		void mcsDisorder2SyncRng(int n);

		virtual double mcsScale() const;

		virtual bool collectH5(splash::DataCollector* data, const std::string& prefix = "");
		virtual bool retrieveH5(splash::DataCollector* data, const std::string& prefix = "");
	};
}
#endif
