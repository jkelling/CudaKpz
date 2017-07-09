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

#ifndef KPZ_KPZ_SIMULATION_CORE_TEMPLATES
#define KPZ_KPZ_SIMULATION_CORE_TEMPLATES

namespace KpzCores
{
	template<class LocalLayout>
	struct SharedBase
	{
		unsigned int system[LocalLayout::BLOCK_DATA];
	};

	__device__ unsigned int normalXor(unsigned int* address, unsigned int val)
	{
		return (*address) ^= val;
	}

	__device__ unsigned int fctAtomicXor(unsigned int* address, unsigned int val)
	{
		return atomicXor(address, val);
	}

	/*! This class uses single-hit doule tiling (DT) DD at block level. Supports DB at device level.
	 */
	template<class Rng, class LocalLayout, class TShared = SharedBase<LocalLayout> >
	class Base : public LocalLayout
	{
		protected:
		int thisIndex, thisShift;
		int xIndex, xShift;
		int yIndex, yShift;
		Rng rng;

		TShared* m_shared;

		public:

		typedef TShared Shared;

			__device__ Base(void* rngData, Shared* shared = 0)
				: rng(rngData), m_shared(shared) {}
		__device__ void setXY(int x, int y) {
			thisIndex = LocalLayout::index(x, y);
			thisShift = LocalLayout::shift(x, y);
			xIndex = LocalLayout::index((x+1), y);
			xShift = LocalLayout::shift((x+1), y);
			yIndex = LocalLayout::index(x, (y+1));
			yShift = LocalLayout::shift(x, (y+1)) | 1; // select slope in y direction
		}

		__device__ int random() {return rng.random();}
		__device__ float randomFloat() {return rng.randomFloat();}

		template<class KpzCore>
		static __device__ void systemLoopDC(KpzCore& core, const int xw, const int yw) {
			const int dispX = (GPU_THREAD_ID_X << (LocalLayout::L_THREAD_CELL_DIM_X_W + LocalLayout::L_BASIC_CELL_DIM_X));
			const int dispY = (GPU_THREAD_ID_Y << (LocalLayout::L_THREAD_CELL_DIM_Y_W + LocalLayout::L_BASIC_CELL_DIM_Y));

			__shared__ int subCell;
			for(int a = 0; a < LocalLayout::UPDATES_PER_THREAD_CELL; ++a) //NOTE number of acive sites is lower
			{
				GPU_LOCAL_SYNC
				int x = core.random();
				if(!(GPU_THREAD_ID_X | GPU_THREAD_ID_Y)) // thread 0 picks reference
				{
					subCell = x>>(LocalLayout::O_WIDTH+LocalLayout::O_WIDTH);
				}
				GPU_LOCAL_SYNC
				const int y = ((x>>LocalLayout::O_WIDTH)&LocalLayout::M_SUB_THREAD_CELL_DIM_Y)
				// (subCell&2) already shift bit by one, thus shift only by (L_SUB_THREAD_CELL_DIM_Y-1)
					+ dispY + ((subCell&2)<<(LocalLayout::L_SUB_THREAD_CELL_DIM_Y-1));
				x = (x&LocalLayout::M_SUB_THREAD_CELL_DIM_X) + dispX + ((subCell&1)<<LocalLayout::L_SUB_THREAD_CELL_DIM_X);

				if((x >= LocalLayout::BLOCK_BORDER_X) || (y >= LocalLayout::BLOCK_BORDER_Y))
				{
#pragma unroll
					for(int s = 0; s < KpzCore::N_RND_SYNC; ++s)
						core.random();
				}
				else
					core.update(x + (xw&(LocalLayout::BASIC_CELL_DIM_X-1)),y + (yw&(LocalLayout::BASIC_CELL_DIM_Y-1)), normalXor);
			}
		}

		__device__ void blockSetupDT(int blocks, int xw, int yw, int block) {}
	};

	/*! This class uses single-hit doule tiling (DT) DD at block level.
	 */
	template<class Rng, class LocalLayout, class TShared = SharedBase<LocalLayout> >
	class BaseDT : public Base<Rng, LocalLayout, TShared>
	{
		typedef Base<Rng, LocalLayout, TShared> MyBase;
		public:

			__device__ BaseDT(void* rngData, typename MyBase::Shared* shared)
				: MyBase(rngData, shared) {}

		template<class KpzCore>
		static __device__ void systemLoopDC(KpzCore& core, const int xw, const int yw) {
			const int dispX = (GPU_THREAD_ID_X << (LocalLayout::L_THREAD_CELL_DIM_X_W + LocalLayout::L_BASIC_CELL_DIM_X));
			const int dispY = (GPU_THREAD_ID_Y << (LocalLayout::L_THREAD_CELL_DIM_Y_W + LocalLayout::L_BASIC_CELL_DIM_Y));

			__shared__ int subCell;
			for(int a = 0; a < LocalLayout::UPDATES_PER_THREAD_CELL; ++a)
			{
				GPU_LOCAL_SYNC
				int x = core.random();
				if(!(GPU_THREAD_ID_X | GPU_THREAD_ID_Y)) // thread 0 picks reference
				{
					subCell = x>>(LocalLayout::O_WIDTH+LocalLayout::O_WIDTH);
				}
				GPU_LOCAL_SYNC
				const int y = ((x>>LocalLayout::O_WIDTH)&LocalLayout::M_SUB_THREAD_CELL_DIM_Y)
					+ dispY + ((subCell&2)<<(LocalLayout::L_SUB_THREAD_CELL_DIM_Y-1));
				x = (x&LocalLayout::M_SUB_THREAD_CELL_DIM_X) + dispX + ((subCell&1)<<LocalLayout::L_SUB_THREAD_CELL_DIM_X);

				core.update(x + (xw&(LocalLayout::BASIC_CELL_DIM_X-1)),y + (yw&(LocalLayout::BASIC_CELL_DIM_Y-1)), normalXor);
			}
		}
	};

#define USING_BASE \
using MyBase::thisIndex; \
using MyBase::thisShift; \
using MyBase::xIndex; \
using MyBase::xShift; \
using MyBase::yIndex; \
using MyBase::yShift; \
using MyBase::setXY; \
using MyBase::m_shared; \
using MyBase::random; \
using MyBase::randomFloat; \
using typename MyBase::Shared; \

	template<class Rng, template<class, class, class> class Base, class LocalLayout >
	class SyncRng
	{
		typedef Base<typename Rng::Device, LocalLayout, SharedBase<LocalLayout> > MyBase;
		void* d_Random;

		public:

		class Device : public MyBase
		{
			public:
			USING_BASE

			static const int N_RND_SYNC = 0; //! Count of random number generated by update() in sync

				__device__ Device(SyncRng& core, Shared* shared)
				: MyBase(core.d_Random, shared) {}

			template<class xor_T>
			inline __device__ void update(int x, int y, const xor_T& myxor) {
				setXY(x,y);

				const int config = ((m_shared->system[thisIndex]>>thisShift)&3)
					| (((m_shared->system[xIndex]>>xShift)&1)<<2) | (((m_shared->system[yIndex]>>yShift)&1)<<3);

				if(config == /*0b1100*/ 12)
				{ // p
					// flip
					myxor(&m_shared->system[thisIndex], 3<<thisShift);
					myxor(&m_shared->system[xIndex], 1<<xShift);
					myxor(&m_shared->system[yIndex], 1<<yShift);
				}
			}

			__device__ void operator() (const int xw, const int yw) {
				MyBase::systemLoopDC(*this, xw, yw);
			}
		};

			SyncRng(Rng& rng, int rndMult)
			: d_Random(rng.getDRandom(LocalLayout::UPDATES_PER_THREAD_CELL*rndMult)) {}
	};

	template<class Rng, template<class, class, class> class Base, class LocalLayout >
	class SyncRngEW : public SyncRng<Rng, Base, LocalLayout>
	{
		typedef SyncRng<Rng, Base, LocalLayout> HostBase;
		typedef typename HostBase::Device MyBase;
		void* d_Random;

		public:

		class Device : public MyBase
		{
			public:
			USING_BASE
			using MyBase::MyBase;

			template<class xor_T>
			inline __device__ void update(int x, int y, const xor_T& myxor) {
				setXY(x,y);

				const int config = ((m_shared->system[thisIndex]>>thisShift)&3)
					| (((m_shared->system[xIndex]>>xShift)&1)<<2) | (((m_shared->system[yIndex]>>yShift)&1)<<3);

				if((config == /*0b1100*/ 12) || (config == /*0b0011*/ 3))
				{ // p = q = 1
					// flip
					myxor(&m_shared->system[thisIndex], 3<<thisShift);
					myxor(&m_shared->system[xIndex], 1<<xShift);
					myxor(&m_shared->system[yIndex], 1<<yShift);
				}
			}

			__device__ void operator() (const int xw, const int yw) {
				MyBase::systemLoopDC(*this, xw, yw);
			}
		};

		using HostBase::HostBase;
	};

	template<class Rng, template<class, class, class> class Base, class LocalLayout>
	class PQSyncRng
	{

		typedef Base<typename Rng::Device, LocalLayout, SharedBase<LocalLayout> > MyBase;
		float p, q;
		void* d_Random;

		public:

		class Device : public MyBase
		{
			float p, q;
			public:
			USING_BASE

			static const int N_RND_SYNC = 1; //! Count of random number generated by update() in sync

				__device__ Device(PQSyncRng& core, Shared* shared)
				: MyBase(core.d_Random, shared), p(core.p), q(core.q) {}

			template<class xor_T>
			__device__ void update(int x, int y, const xor_T& myxor) {
				setXY(x,y);
				const float rnd = randomFloat(); // generate in any case to stay in sync

				const int config = ((m_shared->system[thisIndex]>>thisShift)&3)
					| (((m_shared->system[xIndex]>>xShift)&1)<<2) | (((m_shared->system[yIndex]>>yShift)&1)<<3);

				if(config == /*0b1100*/ 12)
				{ // p
					if(rnd >= p)
						return;
				}
				else if(config == /*0b0011*/ 3)
				{ // q
					if(rnd >= q)
						return;
				}
				else return;

				// flip
				myxor(&m_shared->system[thisIndex], 3<<thisShift);
				myxor(&m_shared->system[xIndex], 1<<xShift);
				myxor(&m_shared->system[yIndex], 1<<yShift);
			}

			__device__ void operator() (const int xw, const int yw) {
				MyBase::systemLoopDC(*this, xw, yw);
			}
		};

			PQSyncRng(Rng& rng, int rndMult, float p, float q)
			: p(p), q(q), d_Random(rng.getDRandom(LocalLayout::UPDATES_PER_THREAD_CELL*(1+Device::N_RND_SYNC)*rndMult)) {}

	};

	template<class Rng, template<class, class, class> class Base, class LocalLayout>
	class Disorder2SyncRng_shared
	{
		struct Shared : public SharedBase<LocalLayout> {
			unsigned int* disorderGlobal;
			unsigned int disorder[LocalLayout::BLOCK_DATA];
		};

		typedef Base<typename Rng::Device, LocalLayout, Shared> MyBase;

		private:
		float p[2], q[2];
		unsigned int* disorderGlobal;
		void* d_Random;

		public:

		class Device : public MyBase
		{
			float p[2], q[2];
			public:
			USING_BASE

			static const int N_RND_SYNC = 1; //! Count of random number generated by update() in sync

				__device__ Device(Disorder2SyncRng_shared& core, Shared* shared)
					: MyBase(core.d_Random, shared) {
					if(GPU_THREAD_ID_X == 0 && GPU_THREAD_ID_Y == 0)
					{
						m_shared->disorderGlobal = core.disorderGlobal;
					}
					p[0] = core.p[0];
					p[1] = core.p[1];
					q[0] = core.q[0];
					q[1] = core.q[1];
				}

			__device__ void blockSetupDT(int blocks, int xw, int yw, int block) {
				copyInDT<LocalLayout>(m_shared->disorder, m_shared->disorderGlobal, blocks, xw>>2, yw>>2, block);
				GPU_LOCAL_SYNC;
			}

			template<class xor_T>
			__device__ void update(int x, int y, const xor_T& myxor) {
#if 1
				setXY(x,y);
				const float rnd = randomFloat(); // generate in any case to stay in sync

				const int config = ((m_shared->system[thisIndex]>>thisShift)&3)
					| (((m_shared->system[xIndex]>>xShift)&1)<<2) | (((m_shared->system[yIndex]>>yShift)&1)<<3);
				/*const int disorderIdx = (y<<L_THREAD_CELL_DIM_X)+x;*/
				/*const int disorderState = (disorder[disorderIdx>>4]>>(disorderIdx&15))&1;*/
				const int disorderState = (m_shared->disorder[thisIndex]>>thisShift)&1;

				if(config == 0b1100)
				{ // p
					if(rnd >= p[disorderState])
						return;
				}
				else if(config == 0b0011)
				{ // q
					if(rnd >= q[disorderState])
						return;
				}
				else return;

				// flip
				myxor(&m_shared->system[thisIndex], 3<<thisShift);
				myxor(&m_shared->system[xIndex], 1<<xShift);
				myxor(&m_shared->system[yIndex], 1<<yShift);
#endif
			}

			__device__ void operator() (const int xw, const int yw) {
				MyBase::systemLoopDC(*this, xw, yw);
			}
		};

			Disorder2SyncRng_shared(Rng& rng, int rndMult, const double* pp, const double* qq, unsigned int* disorderGlobal)
			: disorderGlobal(disorderGlobal)
			, d_Random(rng.getDRandom(LocalLayout::UPDATES_PER_THREAD_CELL*(1+Device::N_RND_SYNC)*rndMult))
		{
			p[0] = pp[0];
			p[1] = pp[1];
			q[0] = qq[0];
			q[1] = qq[1];
		}
	};
		/*! Store disorder for 16 sites in 16 bits of a short int. */
		/*static const int DISORDER_SIZE_W = (LocalLayout::THREAD_CELL_DIM_X_W << LocalLayout::L_THREAD_CELL_DIM_Y_W) >> 1;*/
}

#endif
