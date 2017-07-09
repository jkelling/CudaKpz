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

/*! \file kmcMath.h 
 * \brief Mathematical helper functions.
 * \ingroup DGutils
 */

#ifndef KMC_KMC_MATH_H
#define KMC_KMC_MATH_H

#include <cmath>

/*! Kmc global namespace. */
namespace Kmc
{

/*! \param n integer
 * \return log2 of next highest power of 2 to \p n (or of \p n if it is a power
 * of two)
 */
inline int log2up(unsigned int n)
{
	for(int l = 31; l >= 0; --l, n<<=1)
		if(n&(1<<31))
			if(n<<1)
				return l+1;
			else
				return l;
	return -1;
}

#define FASTLOG2
/*! \param n integer
 * \return integral log2 of \p n
 */
inline int log2(unsigned int n)
{
	#ifndef FASTLOG2
	for(int l = 31; l > 0; --l, n<<=1)
		if(n&(1<<31))
			return l;
	return -1;
	#else //faster version 
	int s = 0;
	for(int a = 16; a > 0; a>>=1)
	{
		if(n >> a)
		{
			n >>= a;
			s +=a;
		}
	}
	return s;
	#endif
}

/*! \param n integer
 * \return log2 of \p n or -1 if \p n is not a power of two
 */
inline int log2Sec(unsigned int n)
{
	if(!n)
		return -1;
	int r = log2(n);
	if(n << (32-r))
		return -1;
	return r;
}

/*! \fn double stdev(double,double,int)
 * \param sum sum of values
 * \param sum2 sum of squares
 * \param n number of values
 * \return standard deviation
 */
inline double stdev(double sum, double sum2, int n)
{
	return sqrt((sum2 - sum*sum/n)/(n-1));
}

}

#endif
