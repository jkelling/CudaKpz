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

#ifndef KMC_DYNAMIC_ORIGIN_DD_H
#define KMC_DYNAMIC_ORIGIN_DD_H

#include "blockDD.h"
#include "dSFMT/dSFMT.h"

namespace Kmc
{
	template<int N>
	class DynamicOriginDD : public BlockDD<N>
	{
		int m_origin[N];
		using BlockDD<N>::m_currentOffset;
		using BlockDD<N>::m_lBlockDim;
		using BlockDD<N>::m_lBlocks;
		using BlockDD<N>::m_mBlocks;
		using BlockDD<N>::m_nBlocks;

		public:
		
			DynamicOriginDD()
		{
			for(int a = 0; a < N; ++a)
				m_origin[a] = 0;
		}

		inline int origin(const int n) const {return m_origin[n];}

		/*! Set specified block current. Calculate the respective offsets.
		 * \param i index of the block
		 */
		inline void setBlock(int i, int blockOffset[N]) const
		{
			for(int a = 0; a < N; ++a)
			{
				blockOffset[a] = ((i&m_mBlocks[a])<<m_lBlockDim[a]) + m_origin[a];
				i >>= m_lBlocks[a];
			}
		}

		/*! Choose random origin.
		 * \param lRestriction origin coordinates have to be multiples of these values.
		 */
		inline void randomOrigin(const int lRestriction[N], dsfmt_t* dsfmt = &dsfmt_global_data)
		{
			for(int a = 0; a < N; ++a)
			{
				const int effBlock = (1<<(m_lBlockDim[a]-lRestriction[a]));
				const int o = int(dsfmt_genrand_close_open(dsfmt)*effBlock);
				m_origin[a] = o<<lRestriction[a];
			}
		}
	};
}

#endif
