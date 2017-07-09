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

#ifndef KPZ_SYSTEMSIZE_H
#define KPZ_SYSTEMSIZE_H

#include "lattice.h"
#include "kpzConst.h"

#include <iostream>

namespace splash
{
	class Dimensions;
}

namespace Kpz
{
	class ThreadLayout;

	struct SystemSizeBase
	{
		// data used on GPU via constmem
		int m_lDimX, m_lDimY;
		int m_dimX, m_dimY;
		int m_lDimXW, m_lDimYW;
		int m_mDimXW, m_mDimYW;

		int m_mBlocksGPU;
		int m_lBlocksX;

		bool set(int lx, int ly);

		inline size_t sizeW() const {return 1<<(m_lDimYW+m_lDimXW);}
	};

	class SystemSize : private SystemSizeBase
	{
		public:

			SystemSize(int lx, int ly) {
				set(lx, ly);
			}

		inline void set(int lx, int ly) {
			SystemSizeBase::set(lx, ly);
		}

		unsigned int* init() const;

		//getters
		inline int lDimX() const {return m_lDimX;}
		inline int lDimY() const {return m_lDimY;}
		inline int dimX() const {return m_dimX;}
		inline int dimY() const {return m_dimY;}
		inline int mDimXW() const {return m_mDimXW;}
		inline int mDimYW() const {return m_mDimYW;}
		inline int lDimXW() const {return m_lDimXW;}
		inline int lDimYW() const {return m_lDimYW;}
		inline int lBlocksX() const {return m_lBlocksX;}
		inline int mBlocksGPU() const {return m_mBlocksGPU;}

		inline int lSize() const {return m_lDimX+m_lDimY;}
		/*! \return size in lattice sites */
		inline unsigned long long size() const {return 1LL<<(lSize());}
		inline unsigned long long sizeW() const {return size()>>4;}
		inline int dimXW() const {return 1<<m_lDimXW;}
		inline int dimYW() const {return 1<<m_lDimYW;}
		inline int blocksGPU() const {return m_mBlocksGPU+1;}

		inline bool validForGPU () const {return m_mBlocksGPU != -1;}
		void adjustDT();

		inline const SystemSizeBase* operator() () const {return this;}

		unsigned int index(int x, int y) const {return (y>>L_BASIC_CELL_DIM_Y<<m_lDimXW) + (x>>L_BASIC_CELL_DIM_X);}
		unsigned int indexW(int xw, int yw) const {return (yw<<m_lDimXW) + xw;}
		static unsigned int shift(int x, int y) {
			return  ((x&3)<<1) // select row by x
					| ((y&3)<<3); // select column by y
		}

		inline LatticePoint latticePoint (int x, int y) const
		{ // reference for latticeencoding
			LatticePoint s;
			s.index = index(x,y);
			s.shift = shift(x,y);
			// lsb of shift selects direction: 0 -> x , 1 -> y
			return s;
		}

		inline Slope slopeX (int x, int y) const
		{
			Slope s;
			s.index = index(x,y);
			s.shift = shift(x,y);
			// lsb of shift selects direction: 0 -> x , 1 -> y
			return s;
		}

		inline Slope slopeY (int x, int y) const
		{
			Slope s;
			s.index = index(x,y);
			s.shift = shift(x,y) | 1; // lsb of shift selects direction: 0 -> x , 1 -> y
			return s;
		}

		inline bool operator==(const SystemSize& b) const {
			return (m_lDimX == b.m_lDimX) & (m_lDimY == b.m_lDimY);
		}

		std::ostream& printConst(std::ostream& o) const;
		void adjust(const ThreadLayout& layout);

		splash::Dimensions splashDimensions() const;

		template<class LocalLayout>
		bool setGPUfromLocalLayout() {
			m_mBlocksGPU = (sizeW()/LocalLayout::BLOCK_SIZE_W) - 1;
			m_lBlocksX = m_lDimXW - LocalLayout::L_BLOCK_DIM_X_W;
			return m_mBlocksGPU==-1;
		}

		unsigned int quadrantIndex_CBBit(int x, int y) const {
			// LSB of x and y identify sublattice, basic cell is 32x1
			return (y<<(m_lDimX-(5+1))) + (x>>(5+1));
		}

		inline Slope slopeX_CBBit (int x, int y) const
		{
			const int sublattice = ((x^y)&1);
			Slope s;
			s.index = quadrantIndex_CBBit(x,y) + ((sizeW()&(0-sublattice)) >> 1);
			s.shift = (x>>1)&31;
			return s;
		}

		inline Slope slopeY_CBBit (int x, int y) const
		{
			const int sublattice = ((x^y)&1);
			const size_t tmp_sizeW = sizeW();
			Slope s;
			s.index = quadrantIndex_CBBit(x,y) + ((((tmp_sizeW<<1)&(0-sublattice)) + tmp_sizeW) >> 2);
			s.shift = (x>>1)&31;
			return s;
		}
	};

	class SystemSizeCPU
	{
		protected:
		int m_mDimX, m_mDimY;

		void cpuSet(const SystemSize& size);

		public:
		
			SystemSizeCPU(const SystemSize& size) {
				cpuSet(size);
			}

		inline void set(const SystemSize& size) {
			cpuSet(size);
		}

		inline int mDimX() const {return m_mDimX;}
		inline int mDimY() const {return m_mDimY;}
	};
	
	std::ostream& operator<<(std::ostream& o, const SystemSize& s);
}

#endif
