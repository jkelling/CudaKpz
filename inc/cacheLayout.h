/***************************************************************************
*   Copyright 2011 - 2012 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef KMC_CACHE_LAYOUT_H
#define KMC_CACHE_LAYOUT_H

#include "kmcMath.h"

#include <fstream>

namespace Kmc
{
	namespace CacheLayoutHelper
	{
		/*! Helper that sums of array components. Array-size templated to allow
		 * loop unrolling by compiler.
		 * \tparam N array size
		 * \param lMax array of int
		 * \return sum of elements of lMax
		 */
		template<int N>
		int totalSize(int lMax[N])
		{
			int sum = lMax[0];
			for(int a = 1; a < N; ++a)
				sum += lMax[a];
			return sum;
		}
	}

	/*! Find a block size and shape that uses at most specified number of bits,
	 * has specified minimum dimensions and is as cubic as possible. Sizes are
	 * given in (simple cubic) lattice sites as log2.
	 * \tparam number of spacial dimensions
	 * \param lMax log2 maximum block size (e.g. system size), also memory to return the result
	 * \param lUseBits log2 desired (maximum) block size
	 * \param lMin log2 minmal allowed block dimensions
	 * \lBitsPerSite log2 number bits per lattice site
	 * \return whether a solution was found
	 */
	template<int N>
	bool getCacheLayout(int lMax[N], const int lUseBits, const int lMin[N], const int lBitsPerSite = 0)
	{
		while((CacheLayoutHelper::totalSize<N>(lMax) + lBitsPerSite) > lUseBits)
		{
			int max = 0, imax = -1;
			for(int a = 0; a < N; ++a)
			{
				if(max < lMax[a] && lMax[a] > lMin[a])
				{
					max = lMax[a];
					imax = a;
				}
			}
			if(imax < 0)
				return false;
			--lMax[imax];
		}
		return true;
	}

	/*! Find a block size and shape that uses at most half the CPU's L1 cache
	 * has specified minimum dimensions and is as cubic as possible. Sizes are
	 * given in (simple cubic) lattice sites as log2.
	 * \tparam number of spacial dimensions
	 * \param lMax log2 maximum block size (e.g. system size), also memory to return the result
	 * \param lMin log2 minmal allowed block dimensions
	 * \lBitsPerSite log2 number bits per lattice site
	 * \return whether a solution was found
	 * \relates Kmc::BlockDD
	 */
	template<int N>
	bool optimizeForL1(int lMax[N], const int lMin[N], const int lBitsPerSite = 0)
	{
		std::ifstream sys("/sys/devices/system/cpu/cpu0/cache/index0/size");
		int bits;
		sys >> bits;
		bits = log2((unsigned)bits) + (10 + 3 - 1); // kb , byte , aim for half the cache size
		return getCacheLayout<N>(lMax, bits, lMin, lBitsPerSite);
	}
}

#endif
