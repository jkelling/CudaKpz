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

#ifndef KPZ_KPZ_SIMULATION_CORE_TEMPLATES
#define KPZ_KPZ_SIMULATION_CORE_TEMPLATES

#include "kpzSimulationPrimitives.cu"

namespace KpzSCACores
{
	template<class KpzCore>
	__device__ void systemLoop(KpzCore& core, unsigned int* system, int parity) {
#pragma unroll
		for(int y = 0; y < 4; ++y)
		{
			core.update(system, (parity^(y&1)), y);
			core.update(system, (parity^(y&1))+2, y);
		}
	}

	template<class Rng>
	class Base
	{
		protected:
		int thisIndex, thisShift;
		int xIndex, xShift;
		int yIndex, yShift;
		unsigned char parity;
		Rng rng;

		__device__ static int index(int x, int y)
		{
			return (y>>Kpz::L_BASIC_CELL_DIM_Y<<1) + (x>>Kpz::L_BASIC_CELL_DIM_X);
		}

		public:

			__device__ Base(void* rngData, unsigned char parity)
				: parity(parity), rng(rngData) {}
		__device__ void setXY(int x, int y) {
			thisIndex = index(x, y);
			thisShift = Kpz::shift(x, y);
			xIndex = index((x+1), y);
			xShift = Kpz::shift((x+1), y);
			yIndex = index(x, (y+1));
			yShift = Kpz::shift(x, (y+1)) | 1; // select slope in y direction
		}

		__device__ int random() {return rng.random();}
		__device__ float randomFloat() {return rng.randomFloat();}

		__device__ inline void initLine(size_t index) {};
		__device__ inline void shiftLine() {};
	};

	template<class Rng>
	class BasePQ : public Base<Rng>
	{
		protected:
		float p, q;

		public:
			__device__ BasePQ(void* rngData, float p, float q, unsigned char parity)
				: Base<Rng>(rngData, parity), p(p), q(q) {}
	};
	
#define USING_BASE \
using Base<typename Rng::Device>::thisIndex; \
using Base<typename Rng::Device>::thisShift; \
using Base<typename Rng::Device>::xIndex; \
using Base<typename Rng::Device>::xShift; \
using Base<typename Rng::Device>::yIndex; \
using Base<typename Rng::Device>::yShift; \
using Base<typename Rng::Device>::setXY; \
using Base<typename Rng::Device>::random; \
using Base<typename Rng::Device>::randomFloat; \
using Base<typename Rng::Device>::parity; \

#define USING_BASEPQ \
using BasePQ<typename Rng::Device>::p; \
using BasePQ<typename Rng::Device>::q; \

	template<class Rng>
	class PQSyncRng
	{
		protected:
		float p, q;
		void* d_Random;

		public:

		class Device : public BasePQ<typename Rng::Device>
		{
			public:
			USING_BASE
			USING_BASEPQ

			static const int N_RND_SYNC = 1; //! Count of random numbers generated per update(x,y) in sync

				__device__ Device(PQSyncRng& core, unsigned char parity)
				: BasePQ<typename Rng::Device>(core.d_Random, core.p, core.q, parity) {}

			__device__ void update(unsigned int* system, int x, int y) {
				setXY(x,y);
				const float rnd = randomFloat(); // generate in any case to stay in sync

				const int config = ((system[thisIndex]>>thisShift)&3)
					| (((system[xIndex]>>xShift)&1)<<2) | (((system[yIndex]>>yShift)&1)<<3);

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
				system[thisIndex] ^= 3<<thisShift;
				system[xIndex] ^= 1<<xShift;
				system[yIndex] ^= 1<<yShift;
			}

			/*! Update system[0] */
			__device__ void update(unsigned int* system) {
				systemLoop(*this, system, parity);
			}
		};

			PQSyncRng(Rng& rng, int rndMult, float p, float q)
			: p(p), q(q), d_Random(rng.getDRandom(Device::N_RND_SYNC*rndMult)) {}
	};

/* Exclude if not needed. Requires CUDA 7.0 */
#ifdef USE_CACHED_RND_OPT
	template<class Rng>
	class PQSyncRngCacheing_5 : public PQSyncRng<Rng>
	{
		typedef PQSyncRng<Rng> MyBase;

		public:

		class Device : public BasePQ<typename Rng::Device>
		{
			int m_rngCache;
			public:
			USING_BASE
			USING_BASEPQ

			static const int N_RND_SYNC = 1; //! Count of random numbers generated per update(x,y) in sync

				__device__ Device(PQSyncRngCacheing_5& core, unsigned char parity)
				: BasePQ<typename Rng::Device>(core.d_Random, core.p, core.q, parity) {}

			__device__ void update(unsigned int* system, int x, int y) {
				setXY(x,y);
				m_rngCache>>=1; // generate in any case to stay in sync

				const int config = ((system[thisIndex]>>thisShift)&3)
					| (((system[xIndex]>>xShift)&1)<<2) | (((system[yIndex]>>yShift)&1)<<3);

				if(config == /*0b1100*/ 12)
				{ // p
					if((m_rngCache&1) >= p)
						return;
				}
				else if(config == /*0b0011*/ 3)
				{ // q
					if((m_rngCache&1) >= q)
						return;
				}
				else return;

				// flip
				system[thisIndex] ^= 3<<thisShift;
				system[xIndex] ^= 1<<xShift;
				system[yIndex] ^= 1<<yShift;
			}

			/*! Update system[0] */
			__device__ void update(unsigned int* system) {
				m_rngCache = random();
				systemLoop(*this, system, parity);
			}
		};

		using MyBase::MyBase;
	};
#endif

__constant__ float DEVICE_CONST_DISORDER2_P[2];
__constant__ float DEVICE_CONST_DISORDER2_Q[2];

	template<class Rng>
	class Disorder2SyncRng
	{
		protected:
		void* d_Random;
		unsigned int* d_Disorder;

		public:

		class Device : public Base<typename Rng::Device>
		{
			unsigned int* d_Disorder;
			unsigned int m_disorder[2];

			public:
			USING_BASE

			static const int N_RND_SYNC = 1; //! Count of random numbers generated per update(x,y) in sync

				__device__ Device(Disorder2SyncRng& core, unsigned char parity)
				: Base<typename Rng::Device>(core.d_Random, parity), d_Disorder(core.d_Disorder) {}

			__device__ void update(unsigned int* system, int x, int y) {
				setXY(x,y);
				const float rnd = randomFloat(); // generate in any case to stay in sync

				const int config = ((system[thisIndex]>>thisShift)&3)
					| (((system[xIndex]>>xShift)&1)<<2) | (((system[yIndex]>>yShift)&1)<<3);

				if(config == /*0b1100*/ 12)
				{ // p
					if(rnd >= DEVICE_CONST_DISORDER2_P[(m_disorder[0]>>thisShift)&1])
						return;
				}
				else if(config == /*0b0011*/ 3)
				{ // q
					if(rnd >= DEVICE_CONST_DISORDER2_Q[(m_disorder[0]>>thisShift)&1])
						return;
				}
				else return;

				// flip
				system[thisIndex] ^= 3<<thisShift;
				system[xIndex] ^= 1<<xShift;
				system[yIndex] ^= 1<<yShift;
			}

			/*! Update system[0] */
			__device__ void update(unsigned int* system) {
				systemLoop(*this, system, parity);
			}

			__device__ inline void initLine(size_t index) {
				m_disorder[0] = d_Disorder[index];
				m_disorder[1] = d_Disorder[index+1];
			}

			__device__ inline void shiftLine() {
				m_disorder[0] = m_disorder[1];
			}
		};

			Disorder2SyncRng(Rng& rng, unsigned int* d_Disorder, int rndMult)
			: d_Random(rng.getDRandom(Device::N_RND_SYNC*rndMult)), d_Disorder(d_Disorder) {}
	};
}
#endif
