/***************************************************************************
*   Copyright 2011 - 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef __OPENCL_VERSION__
#include "kpz.h"

#include "systemSize.h"
#endif

#include "kpzConst.h"

#ifndef __OPENCL_VERSION__
#include <cudaError.h>

using namespace Kpz;
#endif

#include <kmcRandom.h>

#ifndef __OPENCL_VERSION__
__constant__ Kpz::SystemSizeBase KPZ_DEVICE_SYSTEM_SIZE;
	#define SIZE_M_BLOCKS_GPU KPZ_DEVICE_SYSTEM_SIZE.m_mBlocksGPU
	#define SIZE_L_BLOCKS_X KPZ_DEVICE_SYSTEM_SIZE.m_lBlocksX
	#define SIZE_L_DIM_X KPZ_DEVICE_SYSTEM_SIZE.m_lDimX
	#define SIZE_L_DIM_Y KPZ_DEVICE_SYSTEM_SIZE.m_lDimY
	#define SIZE_DIM_X KPZ_DEVICE_SYSTEM_SIZE.m_dimX
	#define SIZE_DIM_Y KPZ_DEVICE_SYSTEM_SIZE.m_dimY
	#define SIZE_L_DIM_X_W KPZ_DEVICE_SYSTEM_SIZE.m_lDimXW
	#define SIZE_L_DIM_Y_W KPZ_DEVICE_SYSTEM_SIZE.m_lDimYW
	#define SIZE_M_DIM_X_W KPZ_DEVICE_SYSTEM_SIZE.m_mDimXW
	#define SIZE_M_DIM_Y_W KPZ_DEVICE_SYSTEM_SIZE.m_mDimYW

#include "cudaDefines.h"

	#define OCL_BASE_LMEM_FKT_CALL
	#define OCL_BASE_LMEM_FKT_HEAD

#include <iostream>
void initStaticDeviceConstDyn(const Kpz::SystemSizeBase* size)
{
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(KPZ_DEVICE_SYSTEM_SIZE, (void*)size, sizeof(Kpz::SystemSizeBase)));
}
#else
	//#define GPU_BLOCK_ID_X (get_global_id(0)>>L_THREADS_X)
	#define GPU_BLOCK_ID_X (get_group_id(0))
	#define GPU_GRID_DIM_X (get_num_groups(0))
	#define GPU_THREAD_ID_X get_local_id(0)
	#define GPU_THREAD_ID_Y get_local_id(1)
	#define GPU_SHARED_MEM __local
	#define GPU_LOCAL_SYNC barrier(CLK_LOCAL_MEM_FENCE);

	// opencl additonal shared mem
	#ifdef KPZ_INNER_DC
		#define BASE_LMEM_SIZE 1
	#else	
		#ifdef WAIT_IF_DEAD
			#define BASE_LMEM_SIZE 4
		#elif defined OPTIMIZED_WAIT
			#define BASE_LMEM_SIZE 2
		#else
			#define BASE_LMEM_SIZE 2
		#endif
	#endif
	#define OCL_BASE_LMEM_FKT_CALL , baseLMem
	#define OCL_BASE_LMEM_FKT_HEAD , __local int* baseLMem

	#define atomicXor atomic_xor
	#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#endif

//#pragma OPENCL EXTENSION cl_amd_printf : enable
#include "stdio.h"

#include "kpzCopy.cu"

template<class KpzCore>
__global__ void kpzKernelTemplate(unsigned int* d_System, const int xw, const int yw, KpzCore core)
{
#if 1 // shared alloc
	// putting this after the first sync, before calculation does not lead to
	// speedup but slowdown in TinyMT case
	GPU_SHARED_MEM typename KpzCore::Device::Shared coreShared;
	typename KpzCore::Device kpzSystem(core, &coreShared);
#endif
	
#if 1 // loop
	int blocks = (GPU_BLOCK_ID_X);
	for(;blocks <= SIZE_M_BLOCKS_GPU; blocks += GPU_GRID_DIM_X)
#endif
	{
#if 1 //copy in
	copyIn<KpzCore::Device>(coreShared.system, d_System, blocks, xw>>2, yw>>2);
	GPU_LOCAL_SYNC
#endif // end copy in

#if 1 // update
	kpzSystem(xw, yw);
	GPU_LOCAL_SYNC
#endif // end update

#if 1 //copy out
	copyOut<KpzCore::Device>(coreShared.system, d_System, blocks, xw>>2, yw>>2);
	GPU_LOCAL_SYNC
#endif //end copy out
	}
}

template<class KpzCore >
__global__ void kpzKernelTemplateDT(unsigned int* d_System, const int xw, const int yw, int block, KpzCore core)
{
#if 1 // shared alloc
	GPU_SHARED_MEM typename KpzCore::Device::Shared coreShared;
	typename KpzCore::Device kpzSystem(core, &coreShared);
#endif
	
#if 1 // loop
	int blocks = (GPU_BLOCK_ID_X);
	for(;blocks <= SIZE_M_BLOCKS_GPU; blocks += GPU_GRID_DIM_X)
#endif
	{
#if 1 //copy in
	copyInDT<KpzCore::Device>(coreShared.system, d_System, blocks, xw>>2, yw>>2, block);
	GPU_LOCAL_SYNC
#endif // end copy in

#if 1 // update
	kpzSystem.blockSetupDT(blocks, xw, yw, block);
	kpzSystem(xw, yw);
	GPU_LOCAL_SYNC
#endif // end update

#if 1 //copy out
	copyOutDT<KpzCore::Device>(coreShared.system, d_System, blocks, xw>>2, yw>>2, block);
	GPU_LOCAL_SYNC
#endif //end copy out
	}
}

// begin Scheduler::callers
#ifndef KPZ_SWITCH_PRAND
#include "kpzRandom.cu"
#else
#include "PrandDevice.cu"
#endif
#include "kpzSimulationPrimitivesTemplates.cu"
#include "scheduler.h"

template<class GpuRng, class LocalLayout>
void Kpz::Scheduler<GpuRng, LocalLayout>::caller_mcsSyncRng(int xw, int yw)
{
	KpzCores::SyncRng<GpuRng, KpzCores::Base, MyLocalLayout> core(m_gpuRng, m_randomNumberBlockMultiplier);
	kpzKernelTemplate <<< m_blocks, MyLocalLayout::threadConfig() >>>
		(d_System, xw, yw, core);
}

template<class GpuRng, class LocalLayout>
void Kpz::Scheduler<GpuRng, LocalLayout>::caller_mcsPQSyncRng(int xw, int yw)
{
	KpzCores::PQSyncRng<GpuRng, KpzCores::Base, MyLocalLayout> core(m_gpuRng, m_randomNumberBlockMultiplier, m_disorderP[0], m_disorderQ[0]);
	kpzKernelTemplate <<< m_blocks, MyLocalLayout::threadConfig() >>>
		(d_System, xw, yw, core);
}

#if KPZ_INNER_DD >= KPZ_INNER_DB && KPZ_INNER_DD < 300 // for future extension
#include "kpzSimulationPrimitivesTemplates_BaseInnerDB.cu"

#if KPZ_INNER_DD == KPZ_INNER_DB
#define Core_BaseDT KpzCores::BaseDTInnerDB
#define Core_BaseDT_KPZ KpzCores::BaseDTInnerNB
#elif KPZ_INNER_DD == KPZ_INNER_DBSH
#define Core_BaseDT KpzCores::BaseDTInnerDBSH
#define Core_BaseDT_KPZ KpzCores::BaseDTInnerNB
#elif KPZ_INNER_DD == KPZ_INNER_DBSHREP
#define Core_BaseDT KpzCores::BaseDTInnerDBSHRep
#define Core_BaseDT_KPZ KpzCores::BaseDTInnerNB
#endif

#elif KPZ_INNER_DD == KPZ_INNER_DT
#define Core_BaseDT KpzCores::BaseDT
#define Core_BaseDT_KPZ KpzCores::BaseDT

#elif KPZ_INNER_DD > KPZ_INNER_DT && KPZ_INNER_DD < KPZ_INNER_DB
#include "kpzSimulationPrimitivesTemplates_BaseInnerDTr.cu"

#define Core_BaseDT KpzCores::BaseDTInnerDTDBSH
#define Core_BaseDT_KPZ KpzCores::BaseDTInnerDTDBSH
#endif

#include "schedulerDT.h"
template<class Scheduler>
void Kpz::SchedulerDTCallers::caller_mcsSyncRng(Scheduler* pthis, int xw, int yw, int block)
{
	typedef typename Scheduler::MyLocalLayout MyLocalLayout;
	KpzCores::SyncRng<decltype(pthis->m_gpuRng), Core_BaseDT_KPZ, MyLocalLayout>
		core(pthis->m_gpuRng, pthis->m_randomNumberBlockMultiplier);
	kpzKernelTemplateDT <<< pthis->m_blocks, MyLocalLayout::threadConfig() >>>
		(pthis->d_System, xw, yw, block, core);
}

template<class Scheduler>
void Kpz::SchedulerDTCallers::caller_mcsPQSyncRng(Scheduler* pthis, int xw, int yw, int block)
{
	typedef typename Scheduler::MyLocalLayout MyLocalLayout;
	if( pthis->m_disorderP[0] == 1. && pthis->m_disorderQ[0] == 1.)
	{
		KpzCores::SyncRngEW<decltype(pthis->m_gpuRng), Core_BaseDT, MyLocalLayout>
			core(pthis->m_gpuRng, pthis->m_randomNumberBlockMultiplier);
		kpzKernelTemplateDT <<< pthis->m_blocks, MyLocalLayout::threadConfig() >>>
			(pthis->d_System, xw, yw, block, core);
	}
	else
	{
		KpzCores::PQSyncRng<decltype(pthis->m_gpuRng), Core_BaseDT, MyLocalLayout>
			core(pthis->m_gpuRng, pthis->m_randomNumberBlockMultiplier, pthis->m_disorderP[0], pthis->m_disorderQ[0]);
		kpzKernelTemplateDT <<< pthis->m_blocks, MyLocalLayout::threadConfig() >>>
			(pthis->d_System, xw, yw, block, core);
	}
}

template<class Scheduler>
void Kpz::SchedulerDTCallers::caller_mcsDisorder2SyncRng(Scheduler* pthis, int xw, int yw, int block)
{
	typedef typename Scheduler::MyLocalLayout MyLocalLayout;
	KpzCores::Disorder2SyncRng_shared<decltype(pthis->m_gpuRng), Core_BaseDT, MyLocalLayout>
		core(pthis->m_gpuRng, pthis->m_randomNumberBlockMultiplier, pthis->m_disorderP, pthis->m_disorderQ, pthis->d_Disorder);
	kpzKernelTemplateDT <<< pthis->m_blocks, MyLocalLayout::threadConfig() >>>
		(pthis->d_System, xw, yw, block, core);
}

namespace Kpz {
	namespace SchedulerDTCallers {
		template void caller_mcsDisorder2SyncRng<SchedulerDT<SLCG64, SchedulerDT_LocalLayoutDis> >
			(SchedulerDT<SLCG64, SchedulerDT_LocalLayoutDis>*, int, int, int);
		template void caller_mcsDisorder2SyncRng<SchedulerDT<TinyMT32, SchedulerDT_LocalLayoutDis> >
			(SchedulerDT<TinyMT32, SchedulerDT_LocalLayoutDis>*, int, int, int);
	}
}

#include "schedulerInstances.h"
#include "schedulerDTInstances.h"
