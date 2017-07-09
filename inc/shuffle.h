/***************************************************************************
*   Copyright 2011 - 2012, 2014 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef KMC_SHUFFLE_H
#define KMC_SHUFFLE_H

#include "helper/kmcTemplateMetaHelper.h"
#include "dSFMT/dSFMT.h"

#include "kmcMath.h"

namespace Kmc
{

/*! \brief Provide permutations of short integer sequences.
 * \ingroup DGutils
 * 
 * \tparam n Number of elements. The actual number of elements will be the nex
 * highest number of two.
 */
template<unsigned int n>
class Shuffle
{
	public:
	static const int L = Log2RoundUp<n>::val;
	static const int N = 1<<L;
	private:
	static const int M = N-1;
	static const int S = 32/L;
	int val[N];

	public:

	/*! \copydetails init() */
	inline Shuffle(dsfmt_t* dsfmt = &dsfmt_global_data) {init(dsfmt);}
	/*! \copydetails init() */
	template<class RandomUnsigned>
	inline Shuffle(RandomUnsigned randomUnsigned) {init(randomUnsigned);}

	/*! Intialize internal array with integer sequence 0...n-1 .*/
	inline void reset() {
		for(int a = 0; a < N; ++a)
			val[a] = a;
	}
	/*! Intialize internal array with integer sequence 0...n-1 and shuffle.*/
	inline void init(dsfmt_t* dsfmt = &dsfmt_global_data) {
		reset();
		shuffle(dsfmt);
	}
	/*! Intialize internal array with integer sequence 0...n-1 and shuffle.*/
	template<class RandomUnsigned>
	inline void init(RandomUnsigned randomUnsigned) {
		reset();
		shuffle(randomUnsigned);
	}

	/*! Permute elements, make use of the method from the key setup of RC4. */
	void shuffle(dsfmt_t* dsfmt = &dsfmt_global_data) {
		for(int a = 0; a < N;)
		{
			unsigned int rnd = (unsigned int)(dsfmt_genrand_close_open(dsfmt)*(4294967296.));
			const int s = a+S;
			for(; a < s && a < N; ++a)
			{
				int i = rnd&M;
				const int v = val[i];
				val[i] = val[a];
				val[a] = v;
				rnd >>= L;
			}
		}
	}
	/*! Permute elements, make use of the method from the key setup of RC4. */
	template<class RandomUnsigned>
	void shuffle(RandomUnsigned randomUnsigned) {
		for(int a = 0; a < N;)
		{
			unsigned int rnd = randomUnsigned();
			const int s = a+S;
			for(; a < s && a < N; ++a)
			{
				int i = rnd&M;
				const int v = val[i];
				val[i] = val[a];
				val[a] = v;
				rnd >>= L;
			}
		}
	}

	/*! \param a index
	 * \return element with index \p a.
	 */
	inline int operator[](int a) const {return val[a];}
};

/*! \brief Provide permutations of short integer sequences.
 * \ingroup DGutils
 * 
 * \details This is a specialization of Shuffle. It allows to the number of
 * elements at runtime.
 */
template<>
class Shuffle<0>
{
	private:
	int L, N, M, S;
	int *val;

	public:
	
	/*! Set the number if elements, round up to next power of two.
	 * \param n number of elements
	 * \return number elements actually allocated
	 */
	int set(int n) {
		const int tL = log2up(n);
		if(tL == L)
		{
			reset();
			return L;
		}
		L = tL;
		N = 1<<L;
		M = N-1;
		S = 32/L;
		delete val;
		val = new int[N];
		reset();
		return L;
	}

	inline Shuffle() : L(-1), val(0) {}
	inline Shuffle(int n, dsfmt_t* dsfmt = &dsfmt_global_data) : L(-1), val(0)
	{
		set(n);
		shuffle(dsfmt);
	}

		~Shuffle() {
		delete[] val;
	}

	/*! Initialize internal array with integer sequence 0...n-1 .*/
	inline void reset() {
		for(int a = 0; a < N; ++a)
			val[a] = a;
	}

	/*! Initialize internal array with integer sequence 0...n-1 and shuffle*/
	inline void init(dsfmt_t* dsfmt = &dsfmt_global_data) {
		reset();
		shuffle(dsfmt);
	}

	/*! Permute elements, make use of the method from the key setup of RC4. */
	void shuffle(dsfmt_t* dsfmt = &dsfmt_global_data) {
		for(int a = 0; a < N;)
		{
			unsigned int rnd = (unsigned int)(dsfmt_genrand_close_open(dsfmt)*(4294967296.));
			const int s = a+S;
			for(; a < s && a < N; ++a)
			{
				int i = rnd&M;
				const int v = val[i];
				val[i] = val[a];
				val[a] = v;
				rnd >>= L;
			}
		}
	}

	/*! \param a index
	 * \return element with index \p a.
	 */
	inline int operator[](int a) const {return val[a];}
	/*! \return log2 number of elements */
	inline int l() const {return L;}
	/*! \return number of elements */
	inline int n() const {return N;}
};

}

#endif
