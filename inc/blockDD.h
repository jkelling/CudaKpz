/***************************************************************************
*   Copyright 2011 - 2012, 2014 - 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef KMC_BLOCK_DD_H
#define KMC_BLOCK_DD_H

#include "shuffle.h"
#include "cacheLayout.h"

namespace Kmc
{
	/*! \brief Handle parameters of block decomposition.
	 * \ingroup DGutils
	 *
	 * \details Handle a block-decomposition of a bit-coded  \c sc lattice (one bit
	 * per site) and random iteration through the blocks. 
	 * \tparam N number of spacial dimensions
	 */
	template<int N>
	class BlockDD
	{
		protected:

		int m_nBlocks, m_blockSites;
		int m_mBlocks[N], m_lBlocks[N], m_lBlockDim[N];
		int m_currentOffset[N], m_currentShift[N];
		Shuffle<0> m_seq;
		int m_index;

		public:

			BlockDD()
		{
			for(int a = 0; a < N; ++a)
				m_currentShift[a] = 0;
		}

		/*! \return number of blocks. */
		inline int nBlocks() const {return m_nBlocks;}
		/*! \return block-size in \c sc lattice sites */
		inline int blockSites() const {return m_blockSites;}
		/*! \param n spacial dimension
		 * \return log2 size of block in direction \p n
		 */
		inline int lBlockDim(const int n) const {return m_lBlockDim[n];}
		/*! \return modulus mask for number of blocks. */
		inline int mBlocks(const int n) const {return m_mBlocks[n];}
		/*! \return log2 number of blocks in dimension n. */
		inline int lBlocks(const int n) const {return m_lBlocks[n];}
		/*! \param n spacial dimension
		 * \return size of block in direction \p n
		 */
		inline int blockDim(const int n) const {return 1<<m_lBlockDim[n];}
		/*! \param n spacial dimension
		 * \return offset of current block in direction \p n
		 */
		inline int blockOffset(const int n) const {return m_currentOffset[n];}
		inline int shift(const int n) const {return m_currentShift[n];}
		/*! \return current index */
		inline int index() const {return m_index;}

		/*! Generate new random sequence of blocks. */
		void shuffle()
		{
			m_index = 0;
			if(m_nBlocks > 1)
			{
				m_seq.shuffle();
				setBlock(m_seq[m_index]);
			}
			else
				setBlock(m_index);
		}
		/*! Generate consecutive sequence of blocks. */
		void resetSeq()
		{
			m_seq.reset();
		}
		/*! Set shift of DD-origin.
		 * \param shift array[N] of shifts, will be truncated to blockDim
		 */
		void setShift(const unsigned int shift[N])
		{
			for(int a = 0; a < N; ++a)
			{
				m_currentShift[a] = shift[a]&((1<<m_lBlockDim[a])-1);
			}
		}

		/*! Make the next block in the sequence current. 
		 * \return whether there is a next block (end of the sequence not
		 * reached)
		 */
		inline bool next()
		{
			++m_index;
			if(m_index >= m_nBlocks)
				return false;
			setBlock(m_seq[m_index]);
			return true;
		}

		/*! Set specified block current. Calculate the respective offsets.
		 * \param i index of the block
		 */
		inline void setBlock(int i)
		{
			for(int a = 0; a < N; ++a)
			{
				m_currentOffset[a] = m_currentShift[a] + ((i&m_mBlocks[a])<<m_lBlockDim[a]);
				i >>= m_lBlocks[a];
			}
		}

		/*! Initialize the decomposition: Determine lateral numbers of blocks, ...
		 * \param layout log2 lateral sizes of blocks in lattice sites for all dimensions
		 * \param size log2 lateral sizes of the system in sites
		 */
		void setLayout(const int layout[N], const int size[N]) {
			m_blockSites = 0;
			m_nBlocks = 0;
			for(int a = 0; a < N; ++a)
			{
				m_blockSites += layout[a];
				m_lBlocks[a] = size[a] - layout[a];
				m_nBlocks += m_lBlocks[a];
				m_mBlocks[a] = (1<<m_lBlocks[a])-1;
				m_lBlockDim[a] = layout[a];
			}
			m_nBlocks = 1 << m_nBlocks;
			m_blockSites = 1 << m_blockSites;
			m_index = 0;
			if(m_nBlocks > 1)
			{
				m_seq.set(m_nBlocks);
				m_seq.init();
				setBlock(m_seq[m_index]);
			}
			else
				setBlock(m_index);
		}

		/*! Choose a layout whre one block uses up around half the L1 cache.
		 * \param size log2 lateral sizes of the system in sites
		 * \param lMin log2 minimum dimensions of blocks
		 */
		void layoutFromL1(const int size[N], const int lMin[N], int lBitsPerSite = 1)
		{
			int lMax[N];
			for(int a = 0; a < N; ++a)
			{
				lMax[a] = size[a];
			}
			Kmc::optimizeForL1<N>(lMax, lMin, lBitsPerSite);
			setLayout(lMax, size);
		}

		static const int N_DIM = N;
	};
}

#endif
