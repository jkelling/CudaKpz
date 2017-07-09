/***************************************************************************
*   Copyright 2014 - 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef KPZ_SCHEDULER_DT_H
#define KPZ_SCHEDULER_DT_H

#include "schedulerBase.h"

#include "kpzConst.h"
#include "LocalLayout.h"

namespace Kpz
{
#ifdef KPZ_FERMI
	typedef LocalLayout<5, 5, 0, 2, 1> SchedulerDT_LocalLayoutDefault;
	typedef LocalLayout<5, 5, 0, 1, 1> SchedulerDT_LocalLayoutDis;
	typedef LocalLayout<5, 5, 2, 1, 1> SchedulerDT_LocalLayoutDis_div2;
	typedef LocalLayout<4, 4, 0, 1, 1> SchedulerDT_LocalLayoutTC1_1T4_4;
	typedef LocalLayout<5, 4, 0, 2, 2> SchedulerDT_LocalLayoutTC2_2;
	typedef LocalLayout<4, 4, 0, 2, 2> SchedulerDT_LocalLayoutTC2_2T4_4;
	typedef LocalLayout<4, 3, 0, 2, 2> SchedulerDT_LocalLayoutTC2_2T4_3;
	typedef LocalLayout<4, 3, 0, 4, 1> SchedulerDT_LocalLayoutTC4_1T4_3;
	typedef LocalLayout<4, 3, 0, 3, 1> SchedulerDT_LocalLayoutTC3_1T4_3;
	typedef LocalLayout<4, 3, 0, 1, 3> SchedulerDT_LocalLayoutTC1_3T4_3;
	typedef LocalLayout<4, 3, 0, 3, 2> SchedulerDT_LocalLayoutTC3_2T4_3;
#else
	typedef LocalLayout<5, 4, 0, 1, 1> SchedulerDT_LocalLayoutDefault;
	typedef LocalLayout<5, 4, 0, 1, 1> SchedulerDT_LocalLayoutDis;
#endif

	namespace SchedulerDTCallers
	{
		template<class Scheduler>
		void caller_mcsSyncRng(Scheduler* pthis, int xw, int yw, int block);
		template<class Scheduler>
		void caller_mcsPQSyncRng(Scheduler* pthis, int xw, int yw, int block);
		template<class Scheduler>
		void caller_mcsDisorder2SyncRng(Scheduler* pthis, int xw, int yw, int block);
	}

	template<class GpuRng, class LocalLayout>
	class SchedulerDT : public SchedulerBase
	{
		typedef SchedulerDT<GpuRng, LocalLayout> Type;
		protected:
		GpuRng m_gpuRng;
		int m_randomNumberBlockMultiplier;
		int calcRandomNumberBlockMultiplier() const {
			int t = m_size.blocksGPU()/m_blocks;
			if(m_size.blocksGPU()%m_blocks > 0)
				++t;
			return t;
		}

		typedef void (*kernelCallerFct) (SchedulerDT<GpuRng, LocalLayout>* pthis, int xw, int yw, int block);
		friend
		void SchedulerDTCallers::caller_mcsSyncRng<>(Type*, int xw, int yw, int block);
		friend
		void SchedulerDTCallers::caller_mcsPQSyncRng<>(Type*, int xw, int yw, int block);
		friend
		void SchedulerDTCallers::caller_mcsDisorder2SyncRng<>(Type*, int xw, int yw, int block);

		template<class A, class B>
		friend struct SchedulerDT_mcsDisorder2SyncRng;

		unsigned int* d_Disorder;

		void doMcs(int n, kernelCallerFct mcs);
		void doMcs(kernelCallerFct mcs);

		void init();
		virtual bool changedSizeEvent();

		public:
		typedef LocalLayoutProxy<LocalLayout, LocalLayoutDT> MyLocalLayout;

			SchedulerDT(int lx, int ly, int device = -1);
			SchedulerDT(const SystemSize& size, int device = -1);
			SchedulerDT(SchedulerService &&other, int device = -1);
			~SchedulerDT();

		void randomize() {m_gpuRng.randomize();}

		void mcs(int n) {
			doMcs(n, &SchedulerDTCallers::caller_mcsSyncRng<Type>);
		}
		void mcsSyncRng(int n) {
			doMcs(n, &SchedulerDTCallers::caller_mcsSyncRng<Type>);
		}

		void mcsPQSyncRng(int n) {
			doMcs(n, &SchedulerDTCallers::caller_mcsPQSyncRng<Type>);
		}
		void mcsDisorder2SyncRng(int n);

		void releaseDeviceDisorder();

		virtual bool collectH5(splash::DataCollector* data, const std::string& prefix = "");
		virtual bool retrieveH5(splash::DataCollector* data, const std::string& prefix = "");
	};
}
#endif
