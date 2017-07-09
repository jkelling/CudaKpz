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

#ifndef KMC_DOUBLE_TILING_DD_H
#define KMC_DOUBLE_TILING_DD_H

#include "blockDD.h"

#include <iostream>

namespace Kmc
{
	/*!	\brief View to handle logical double tiling (DT) domain decompostion
	 *
	 * Size of sub-blocks is given by properties of BlockDD base. The number of
	 * blocks in BlockDD must be even for DT to work.
	 */
	template<int N>
	class DoubleTilingDD : public BlockDD<N>
	{
		int m_domainSet;
		using BlockDD<N>::m_lBlockDim;
		using BlockDD<N>::m_lBlocks;
		using BlockDD<N>::m_mBlocks;
		using BlockDD<N>::m_currentOffset;
		using BlockDD<N>::m_currentShift;

		public:
		
			DoubleTilingDD() : m_domainSet(0)
		{
		}

		inline int domainSet() const {return m_domainSet;}
		inline void setDomainSet(int d) {m_domainSet = d;}
		inline int nBlocks() const {return BlockDD<N>::nBlocks()>>N;}

		/*! Set specified block current. Calculate the respective offsets.
		 * \param i index of the block
		 */
		inline void setBlock(int i)
		{
			for(int a = 0; a < N; ++a)
			{
				m_currentOffset[a] = m_currentShift[a] + ((((i<<=1)&m_mBlocks[a]) + ((m_domainSet>>a)&1))<<(m_lBlockDim[a]));
				i >>= m_lBlocks[a];
			}
		}
		/*! Set specified block current. Calculate the respective offsets.
		 * \param i index of the block
		 * \param blockOffset arrax[N] to return block offsets
		 */
		inline void setBlock(int i, int blockOffset[N]) const
		{
			for(int a = 0; a < N; ++a)
			{
				blockOffset[a] = m_currentShift[a] + ((((i<<=1)&m_mBlocks[a]) + ((m_domainSet>>a)&1))<<(m_lBlockDim[a]));
				i >>= m_lBlocks[a];
			}
		}

		static const int N_DOMAIN_SET = 1<<N;
	};
}

#endif
