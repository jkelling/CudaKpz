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

namespace KpzSCABitCores
{
	template<class Rng>
	class Base
	{
		protected:
		Rng rng;

		public:

			__device__ Base(void* rngData)
				: rng(rngData) {}

		__device__ int random() {return rng.random();}
		__device__ float randomFloat() {return rng.randomFloat();}
	};

	template<class Rng>
	class Empty
	{
		void* d_Random;

		public:

		class Device : public Base<typename Rng::Device>
		{
			typedef Base<typename Rng::Device> MyBase;
			public:
			static const int N_RND_SYNC = 1; //! Count of random numbers generated per update(x,y) in sync

				__device__ Device(Empty& core)
				: Base<typename Rng::Device>(core.d_Random) {}

			using MyBase::random;

			class UpdateMask
			{
				protected:
				unsigned int m_mask;

				__device__ inline UpdateMask(unsigned int mask) : m_mask(mask) {}
				public:

				__device__ inline UpdateMask(Device& dev) : m_mask(dev.random()) {}

				__device__ inline void genMask
					(const unsigned int thisDataX, const unsigned int thisDataY, const unsigned int xDataX, const unsigned int yDataY) {
					// both x and y of this unset
					// ... and both xx and yy set
					m_mask &= ~(thisDataX | thisDataY) & xDataX & yDataY;
				}
				__device__ unsigned int mask() const {return m_mask;}
			};
		};

			Empty(Rng& rng, int rndMult)
			: d_Random(rng.getDRandom(Device::N_RND_SYNC*rndMult)) {}
	};

	/*template<class Rng, template<class Dev> unsigned int (*GenP)(Dev&)>*/
	template<class Rng, unsigned int (*genP)(Base<typename Rng::Device>&)>
	class EmptyGenP : public Empty<Rng>
	{
		public:

		typedef typename Empty<Rng>::Device DevBase;
		class Device : public DevBase
		{
			public:
			using DevBase::DevBase;
			using DevBase::random;

			typedef typename DevBase::UpdateMask UpBase;
			class UpdateMask : public UpBase
			{
				protected:
				using UpBase::m_mask;
				using UpBase::UpBase;

				public:

				__device__ inline UpdateMask(Device& dev) : UpBase(genP(dev)) {}

			};
		};

		using Empty<Rng>::Empty;
	};

	/*template<class Rng, template<class Dev> unsigned int (*GenPQ)(Dev&)>*/
	template<class Rng, unsigned int (*genP)(Base<typename Rng::Device>&)>
	class EmptyGenPQ : public EmptyGenP<Rng, genP>
	{
		public:

		typedef typename EmptyGenP<Rng, genP>::Device DevBase;
		class Device : public DevBase
		{
			public:
			using DevBase::DevBase;
			using DevBase::random;

			typedef typename DevBase::UpdateMask UpBase;
			class UpdateMask : public UpBase
			{
				protected:
				unsigned int m_maskQ;
				using UpBase::m_mask;
				using UpBase::UpBase;

				public:
				static const int N_RND_SYNC = 2; //! Count of random numbers generated per update(x,y) in sync

				__device__ inline UpdateMask(Device& dev) : UpBase(dev), m_maskQ(genP(dev)) {}

				__device__ inline void genMask
					(const unsigned int thisDataX, const unsigned int thisDataY, const unsigned int xDataX, const unsigned int yDataY) {
					UpBase::genMask(thisDataX, thisDataY, xDataX, yDataY);
					// both x and y of this set
					// ... and both xx and yy unset
					m_mask ^= m_maskQ & ~(xDataX | yDataY) & thisDataX & thisDataY;
				}
			};
		};

		using EmptyGenP<Rng, genP>::EmptyGenP;
	};


	template<class Rng>
	class EmptyPQ : public Empty<Rng>
	{
		public:

		typedef typename Empty<Rng>::Device DevBase;
		class Device : public DevBase
		{
			public:
			static const int N_RND_SYNC = 2; //! Count of random numbers generated per update(x,y) in sync

			using DevBase::DevBase;
			using DevBase::random;

			typedef typename DevBase::UpdateMask UpBase;
			class UpdateMask : public UpBase
			{
				protected:
				unsigned int m_maskQ;
				using UpBase::m_mask;

				__device__ inline UpdateMask(unsigned int mask) : UpBase(mask), m_maskQ(mask) {}
				public:

				__device__ inline UpdateMask(Device& dev) : UpBase(dev), m_maskQ(dev.random()) {}

				__device__ inline void genMask
					(const unsigned int thisDataX, const unsigned int thisDataY, const unsigned int xDataX, const unsigned int yDataY) {
					UpBase::genMask(thisDataX, thisDataY, xDataX, yDataY);
					// both x and y of this set
					// ... and both xx and yy unset
					m_mask ^= m_maskQ & ~(xDataX | yDataY) & thisDataX & thisDataY;
				}
			};
		};

			EmptyPQ(Rng& rng, int rndMult)
			: Empty<Rng>(rng, rndMult) {}
	};

	__shared__ float sharedPQ[4];

	template<class Rng>
	class EmptyAnyP : public Empty<Rng>
	{
		float m_p;

		public:

		typedef typename Empty<Rng>::Device DevBase;
		class Device : public DevBase
		{
			public:

				__device__ Device(EmptyAnyP& core)
				: DevBase(core) {
					if(threadIdx.x == 0)
						sharedPQ[0] = core.m_p;
				}

			using DevBase::randomFloat;

			typedef typename DevBase::UpdateMask UpBase;
			class UpdateMask : public UpBase
			{
				protected:
				using UpBase::m_mask;

				public:

				__device__ inline UpdateMask(Device& dev) : UpBase(0) {
					for(int a = 0; a < 32; ++a)
						m_mask |= (dev.randomFloat()<sharedPQ[0])<<a;
				}
			};
		};

			EmptyAnyP(Rng& rng, int rndMult, float p)
			: Empty<Rng>(rng, rndMult), m_p(p) {}
	};

	template<class Rng>
	class EmptyAnyPQ : public EmptyPQ<Rng>
	{
		float m_pq[2];

		public:

		typedef typename EmptyPQ<Rng>::Device DevBase;
		class Device : public DevBase
		{
			public:

				__device__ Device(EmptyAnyPQ& core)
				: DevBase(core) {
					if(threadIdx.x < 2)
						sharedPQ[threadIdx.x] = core.m_pq[threadIdx.x];
				}

			using DevBase::randomFloat;

			typedef typename DevBase::UpdateMask UpBase;
			class UpdateMask : public UpBase
			{
				protected:
				using UpBase::m_mask;
				using UpBase::m_maskQ;

				public:

				__device__ inline UpdateMask(Device& dev) : UpBase(0) {
					for(int a = 0; a < 32; ++a)
					{
						m_mask |= (dev.randomFloat()<sharedPQ[0])<<a;
						m_maskQ |= (dev.randomFloat()<sharedPQ[1])<<a;
					}
				}
			};
		};

			EmptyAnyPQ(Rng& rng, int rndMult, float p, float q)
			: EmptyPQ<Rng>(rng, rndMult) {
			m_pq[0] = p;
			m_pq[1] = q;
		}
	};
}
