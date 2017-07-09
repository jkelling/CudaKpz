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

namespace KpzCores
{
	/*! This class uses dead border (DB) DD with delayed borders at block level.
	 */
	template<class Rng, class LocalLayout, class TShared = SharedBase<LocalLayout> >
	class BaseDTInnerDB : public BaseDT<Rng, LocalLayout, TShared>
	{
		static_assert(LocalLayout::L_THREAD_CELL_DIM_X <= 8, "Thread cells are too large in x direction (max is 256).");
		static_assert(LocalLayout::L_THREAD_CELL_DIM_Y <= 8, "Thread cells are too large in y direction (max is 256).");

		typedef BaseDT<Rng, LocalLayout, TShared> MyBase;
		public:

			__device__ BaseDTInnerDB(void* rngData, typename MyBase::Shared* shared)
				: MyBase(rngData, shared) {}

		template<class KpzCore>
		static __device__ void systemLoopDC(KpzCore& core, const int xw, const int yw) {
			const int dispX = (GPU_THREAD_ID_X << (LocalLayout::L_THREAD_CELL_DIM_X_W + LocalLayout::L_BASIC_CELL_DIM_X))
				+ (xw&(LocalLayout::BASIC_CELL_DIM_X-1));
			const int dispY = (GPU_THREAD_ID_Y << (LocalLayout::L_THREAD_CELL_DIM_Y_W + LocalLayout::L_BASIC_CELL_DIM_Y))
				+ (yw&(LocalLayout::BASIC_CELL_DIM_Y-1));

			for(int a = 0; a < LocalLayout::UPDATES_PER_THREAD_CELL; ++a)
			{
				int x = core.random();
				int y = ((x>>LocalLayout::O_WIDTH)&LocalLayout::M_THREAD_CELL_DIM_Y);
				x = (x&LocalLayout::M_THREAD_CELL_DIM_X);
				const bool border = (y == LocalLayout::M_THREAD_CELL_DIM_Y | x == LocalLayout::M_THREAD_CELL_DIM_X);
				const bool borderCorner = (y == LocalLayout::M_THREAD_CELL_DIM_Y & x == LocalLayout::M_THREAD_CELL_DIM_X);
				x += dispX;
				y += dispY;
				x &= LocalLayout::M_BLOCK_DIM_X;
				y &= LocalLayout::M_BLOCK_DIM_Y;
						
				GPU_LOCAL_SYNC
				if(!border)
					core.update(x, y, fctAtomicXor);
				GPU_LOCAL_SYNC
				if(border != borderCorner)
					core.update(x, y, fctAtomicXor);
				GPU_LOCAL_SYNC
				if(borderCorner)
					core.update(x, y, fctAtomicXor);
			}
		}
	};

	/*! This class uses block DD with no borders at block level. This is
	 * equivalent to any type of DB if p or q is 0. It will produce
	 * inconsitencies if p!=0 and q!=0.
	 */
	template<class Rng, class LocalLayout, class TShared = SharedBase<LocalLayout> >
	class BaseDTInnerNB : public BaseDTInnerDB<Rng, LocalLayout, TShared>
	{
		typedef BaseDTInnerDB<Rng, LocalLayout, TShared> MyBase;
		public:

			__device__ BaseDTInnerNB(void* rngData, typename MyBase::Shared* shared)
				: MyBase(rngData, shared) {}

		template<class KpzCore>
		static __device__ void systemLoopDC(KpzCore& core, const int xw, const int yw) {
			const int dispX = (GPU_THREAD_ID_X << (LocalLayout::L_THREAD_CELL_DIM_X_W + LocalLayout::L_BASIC_CELL_DIM_X))
				+ (xw&(LocalLayout::BASIC_CELL_DIM_X-1));
			const int dispY = (GPU_THREAD_ID_Y << (LocalLayout::L_THREAD_CELL_DIM_Y_W + LocalLayout::L_BASIC_CELL_DIM_Y))
				+ (yw&(LocalLayout::BASIC_CELL_DIM_Y-1));

			for(int a = 0; a < LocalLayout::UPDATES_PER_THREAD_CELL; ++a)
			{
				int x = core.random();
				int y = ((x>>LocalLayout::O_WIDTH)&LocalLayout::M_THREAD_CELL_DIM_Y) + dispY;
				x = (x&LocalLayout::M_THREAD_CELL_DIM_X) + dispX;
				x &= LocalLayout::M_BLOCK_DIM_X;
				y &= LocalLayout::M_BLOCK_DIM_Y;

				GPU_LOCAL_SYNC
				core.update(x, y, fctAtomicXor);
			}
		}
	};

	/*! This class uses single-hit dead border (DBSH) DD with delayed borders at block level.
	 */
	template<class Rng, class LocalLayout, class TShared = SharedBase<LocalLayout> >
	class BaseDTInnerDBSH : public BaseDTInnerDB<Rng, LocalLayout, TShared>
	{
		protected:
		typedef BaseDTInnerDB<Rng, LocalLayout, TShared> MyBase;

		union Origin {
			short value;
			struct {
				char x, y;
			};
		};

		public:

			__device__ BaseDTInnerDBSH(void* rngData, typename MyBase::Shared* shared)
				: MyBase(rngData, shared) {}

		template<class KpzCore>
		static __device__ void systemLoopDC(KpzCore& core, const int xw, const int yw) {
			const int dispX = (GPU_THREAD_ID_X << (LocalLayout::L_THREAD_CELL_DIM_X_W + LocalLayout::L_BASIC_CELL_DIM_X))
				+ (xw&(LocalLayout::BASIC_CELL_DIM_X-1));
			const int dispY = (GPU_THREAD_ID_Y << (LocalLayout::L_THREAD_CELL_DIM_Y_W + LocalLayout::L_BASIC_CELL_DIM_Y))
				+ (yw&(LocalLayout::BASIC_CELL_DIM_Y-1));

			__shared__ Origin origin;
			for(int a = 0; a < LocalLayout::UPDATES_PER_THREAD_CELL; ++a)
			{
				int x = core.random();
				if(!(GPU_THREAD_ID_X | GPU_THREAD_ID_Y)) // thread 0 picks reference
				{
					origin.value = x>>16;
				}
				GPU_LOCAL_SYNC
				int y = ((x>>LocalLayout::O_WIDTH)&LocalLayout::M_THREAD_CELL_DIM_Y);
				x = (x&LocalLayout::M_THREAD_CELL_DIM_X);
				const bool border = (y == LocalLayout::M_THREAD_CELL_DIM_Y | x == LocalLayout::M_THREAD_CELL_DIM_X);
				const bool borderCorner = (y == LocalLayout::M_THREAD_CELL_DIM_Y & x == LocalLayout::M_THREAD_CELL_DIM_X);
				x += dispX + (origin.x&LocalLayout::M_THREAD_CELL_DIM_X);
				y += dispY + (origin.y&LocalLayout::M_THREAD_CELL_DIM_Y);
				x &= LocalLayout::M_BLOCK_DIM_X;
				y &= LocalLayout::M_BLOCK_DIM_Y;

				if(!border)
					core.update(x, y, fctAtomicXor);
				GPU_LOCAL_SYNC
				if(border != borderCorner)
					core.update(x, y, fctAtomicXor);
				GPU_LOCAL_SYNC
				if(borderCorner)
					core.update(x, y, fctAtomicXor);
			}
		}
	};

	/*! This class uses single-hit dead border (DB) DD at block level.
	 * 
	 * The preformance is about 10% increased over BaseDTInnerDB replaceing all
	 * border hitting updates until after the main update loop. At this point
	 * border hitting updates are carried out in sequence, this could be
	 * further optimized adding more iterations of replaced bordr hitting
	 * updates.
	 * \todo Add adjusting number of delay iterations if this is to be used with
	 * small TC.
	 */
	template<class Rng, class LocalLayout, class TShared = SharedBase<LocalLayout> >
	class BaseDTInnerDBSHRep : public BaseDTInnerDBSH<Rng, LocalLayout, TShared>
	{
		protected:
		typedef BaseDTInnerDBSH<Rng, LocalLayout, TShared> MyBase;
		using typename MyBase::Origin;

		public:
			__device__ BaseDTInnerDBSHRep(void* rngData, typename MyBase::Shared* shared)
				: MyBase(rngData, shared) {}

		template<class KpzCore>
		static __device__ void systemLoopDC(KpzCore& core, const int xw, const int yw) {
			const int dispX = (GPU_THREAD_ID_X << (LocalLayout::L_THREAD_CELL_DIM_X_W + LocalLayout::L_BASIC_CELL_DIM_X))
				+ (xw&(LocalLayout::BASIC_CELL_DIM_X-1));
			const int dispY = (GPU_THREAD_ID_Y << (LocalLayout::L_THREAD_CELL_DIM_Y_W + LocalLayout::L_BASIC_CELL_DIM_Y))
				+ (yw&(LocalLayout::BASIC_CELL_DIM_Y-1));

			__shared__ Origin origin;
			unsigned int borderHits = 0;
			for(int a = 0; a < LocalLayout::UPDATES_PER_THREAD_CELL; ++a)
			{
				int x = core.random();
				if(!(GPU_THREAD_ID_X | GPU_THREAD_ID_Y)) // thread 0 picks reference
				{
					origin.value = x>>16;
				}
				GPU_LOCAL_SYNC
				int y = ((x>>LocalLayout::O_WIDTH)&LocalLayout::M_THREAD_CELL_DIM_Y);
				x = (x&LocalLayout::M_THREAD_CELL_DIM_X);
				const bool border = (y == LocalLayout::M_THREAD_CELL_DIM_Y | x == LocalLayout::M_THREAD_CELL_DIM_X);
				if(border)
					++borderHits;
				else
				{
					x += dispX + (origin.x&LocalLayout::M_THREAD_CELL_DIM_X);
					y += dispY + (origin.y&LocalLayout::M_THREAD_CELL_DIM_Y);
					x &= LocalLayout::M_BLOCK_DIM_X;
					y &= LocalLayout::M_BLOCK_DIM_Y;
						
					core.update(x, y, fctAtomicXor);
				}
				GPU_LOCAL_SYNC
			}
			// warp-reduce borderHits
			__shared__ unsigned int sBorderHits;
			if(!(GPU_THREAD_ID_X | GPU_THREAD_ID_Y))
				sBorderHits = borderHits;
			int id = (GPU_THREAD_ID_Y << LocalLayout::L_THREADS_X) + GPU_THREAD_ID_X;
#if __CUDA_ARCH__ >= 300
			for(int a = 16; a >= 1; a>>=1)
				borderHits += __shfl_xor(borderHits, a, 32);
			GPU_LOCAL_SYNC
			if(((id & 31) == 0) && id > 0)
#endif
			// Just have all threads hammer the memory with the atomic if we
			// cannot do warp shuffle. This is ineffiecient, but at least it
			// will compile for Fermi
				atomicAdd(&sBorderHits, borderHits);
			GPU_LOCAL_SYNC
			for(int a = 0 ; a < sBorderHits; a += LocalLayout::THREADS)
			{
				const bool mask = (a+id) < sBorderHits;
				int x;
				if(mask)
				{
					x = core.random();
					if(!(GPU_THREAD_ID_X | GPU_THREAD_ID_Y)) // thread 0 picks reference
					{
						origin.value = x>>16;
					}
				}
				GPU_LOCAL_SYNC
				int y = ((x>>LocalLayout::O_WIDTH)&LocalLayout::M_THREAD_CELL_DIM_Y);
				x = (x&LocalLayout::M_THREAD_CELL_DIM_X);
				const bool border = (y == LocalLayout::M_THREAD_CELL_DIM_Y | x == LocalLayout::M_THREAD_CELL_DIM_X);
				const bool borderCorner = (y == LocalLayout::M_THREAD_CELL_DIM_Y & x == LocalLayout::M_THREAD_CELL_DIM_X);
				if(mask)
				{
					x += dispX + (origin.x&LocalLayout::M_THREAD_CELL_DIM_X);
					y += dispY + (origin.y&LocalLayout::M_THREAD_CELL_DIM_Y);
					x &= LocalLayout::M_BLOCK_DIM_X;
					y &= LocalLayout::M_BLOCK_DIM_Y;
							
					if(!border)
						core.update(x, y, fctAtomicXor);
				}
				GPU_LOCAL_SYNC
				if(border & mask & !borderCorner)
					core.update(x, y, fctAtomicXor);
				GPU_LOCAL_SYNC
				if(mask & borderCorner)
					core.update(x, y, fctAtomicXor);
			}
		}
	};
}
