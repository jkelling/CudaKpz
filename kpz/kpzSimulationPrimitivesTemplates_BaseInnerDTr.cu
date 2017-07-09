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

#include "kpzSimulationPrimitivesTemplates.cu"
#include "kpzSimulationPrimitivesTemplates_BaseInnerDB.cu"

namespace KpzCores
{

	/*! This class uses single-hit doule tiling (DT) DD with random origin at block level.
	 */
	template<class Rng, class LocalLayout, class TShared = SharedBase<LocalLayout> >
	class BaseDTInnerDTDBSH : public BaseDTInnerDBSH<Rng, LocalLayout, TShared>
	{
		static_assert(!(LocalLayout::UPDATES_PER_THREAD_CELL % 32), "UPDATES_PER_THREAD_CELL must be multiple of warpsize (32).");
		static_assert(LocalLayout::THREADS >= 32, "Need to have a minimum of 32 threads per block (one warp).");

		typedef BaseDTInnerDBSH<Rng, LocalLayout, TShared> MyBase;
		using typename MyBase::Origin;
		
		protected:

		union OriginSub {
			unsigned int value;
			struct {
				short subCell;
				Origin origin;
			};
		};
		static const unsigned int ORIGIN_MASK = 0xFFFF
			| (LocalLayout::M_THREAD_CELL_DIM_X<<16) | (LocalLayout::M_THREAD_CELL_DIM_Y<<(LocalLayout::O_WIDTH+16));

		public:

			__device__ BaseDTInnerDTDBSH(void* rngData, typename MyBase::Shared* shared)
				: MyBase(rngData, shared) {}

		template<class KpzCore>
		static __device__ void systemLoopDC(KpzCore& core, const int xw, const int yw) {
			const int dispX = (GPU_THREAD_ID_X << (LocalLayout::L_THREAD_CELL_DIM_X_W + LocalLayout::L_BASIC_CELL_DIM_X));
			const int dispY = (GPU_THREAD_ID_Y << (LocalLayout::L_THREAD_CELL_DIM_Y_W + LocalLayout::L_BASIC_CELL_DIM_Y));

			__shared__ OriginSub o[32];
			for(int a = 0; a < LocalLayout::UPDATES_PER_THREAD_CELL; ++a)
			{
				if(! (a&31) )
				{
					const int tid = (GPU_THREAD_ID_Y << LocalLayout::L_THREADS_X) + GPU_THREAD_ID_X;
					if(tid < 32)
						o[tid].value = core.random()&ORIGIN_MASK;
				}
				GPU_LOCAL_SYNC
				int x = core.random();
				int y = ((x>>LocalLayout::O_WIDTH)&LocalLayout::M_SUB_THREAD_CELL_DIM_Y) + dispY;
				x = (x&LocalLayout::M_SUB_THREAD_CELL_DIM_X) + dispX;

				x += ((o[a&31].subCell&1)<<LocalLayout::L_SUB_THREAD_CELL_DIM_X) + o[a&31].origin.x;
				y += ((o[a&31].subCell&2)<<(LocalLayout::L_SUB_THREAD_CELL_DIM_Y-1)) + o[a&31].origin.y;
				x &= LocalLayout::M_BLOCK_DIM_X;
				y &= LocalLayout::M_BLOCK_DIM_Y;

				core.update(x + (xw&(LocalLayout::BASIC_CELL_DIM_X-1)),y + (yw&(LocalLayout::BASIC_CELL_DIM_Y-1)), normalXor);
			}
		}
	};
}
