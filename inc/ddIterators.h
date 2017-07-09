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

#ifndef KMC_DYNAMIC_DD_ITERATORS_H
#define KMC_DYNAMIC_DD_ITERATORS_H

namespace Kmc
{
	template<class DD>
	class DDIterator
	{
		protected:
		DD& m_dd;
		unsigned int m_index;
		int m_currentOffset[DD::N_DIM];

		public:
			DDIterator(DD& dd, int index = 0) : m_dd(dd) {setIndex(index);}
			DDIterator(DDIterator& other) : m_dd(other.m_dd), m_index(other.m_index)
		{
			for(int a = 0; a < DD::N_DIM; ++a)
				m_currentOffset[a] = other.m_currentOffset[a];
		}

		/*! Set index to i and make block current. Do nothing if i is out of range.
		 * \param i block index
		 * \return success
		 */
		inline bool setIndex(unsigned int i)
		{
			if(i >= m_dd.nBlocks())
				return false;
			m_index = i;
			m_dd.setBlock(m_index, m_currentOffset);
			return true;
		}

		/*! \param n spacial dimension
		 * \return offset of current block in direction \p n
		 */
		inline int blockOffset(const int n) const {return m_currentOffset[n];}
		inline int index() const {return m_index;}

		/*! Make the next block in the sequence current. 
		 * \return whether there is a next block (end of the sequence not
		 * reached)
		 */
		inline bool next(int increment = 1)
		{
			return setIndex(m_index+increment);
		}
	};
}

#endif
