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

/*! \file kmcRandom.h 
 * \brief LCRNG parameters.
 *
 * This file contains parameters for both an 32-bit LCRNG and the strided
 * 64-bit LCG used in GPU code.
 */

#ifndef KMC_RANDOM_H
#define KMC_RANDOM_H

#ifndef KMC_NO_RNG_SHUFFLING
#ifndef KMC_L_LCG_N
#warning "KMC_L_LCG_N not set, setting it to 4"
#define KMC_L_LCG_N 4
#endif
#endif

const int KMC_L_LCG_NATIVE_PERIOD = 30; //save assumption, should in fact be 2^31-1
const int KMC_LCG_A = 16807;
const int KMC_LCG_M = 2147483647;
const int KMC_LCG_Q = 127773;
const int KMC_LCG_R = 2836;
#ifndef KMC_NO_RNG_SHUFFLING
const int KMC_LCG_N = (1<<KMC_L_LCG_N);
//for this to work KMC_LCG_N has to be a power of 2
const int KMC_LCG_N_M = (KMC_LCG_N-1);
const int KMC_LCG_DWORDS = KMC_LCG_N + 2; // plus two for k and i
#else
const int KMC_LCG_DWORDS = 1; // only i
#endif

const float KMC_LCG_RAND_SUP = KMC_LCG_M - 1.;
const float KMC_LCG_RAND_SUP_RED = 1<<23; 
const int KMC_M_LCG_RAND_RED = (int)KMC_LCG_RAND_SUP_RED-1;

// 64bit strided LCG
#ifndef __OPENCL_VERSION__
const unsigned long long KMC_SLCG_A = 2862933555777941757ull;
const unsigned long long KMC_SLCG_C = 1442695040888963407ull;
#else
#define KMC_SLCG_A = 2862933555777941757lu;
#define KMC_SLCG_C = 1442695040888963407lu;
#endif

#ifndef __OPENCL_VERSION__
// the following two functions do not work inline for static const class members
// and for CUDA:
/*! Compute modified parameter \c A for the 64-bit LCG to skip ahead in the
 * sequence.
 * \param n number of generators
 * \return modified parameter \c A
 */
inline static unsigned long long SLCGskipA (int n) {
	long long a = KMC_SLCG_A;
	for(; n>1; --n) a *= KMC_SLCG_A;
	return a;
}
/*! Compute modified parameter \c C for the 64-bit LCG to skip ahead in the
 * sequence.
 * \param n number of generators
 * \return modified parameter \c C
 */
inline static unsigned long long SLCGskipC (int n) {
	long long a = 1, c = 0;
	for(; n>0; --n) {
		a *= KMC_SLCG_A;
		c += a*KMC_SLCG_C;
	}
	return c;
}

// alternative: template metaprogramming
#include "helper/kmcTemplateMetaHelper.h"
/*! \brief Template with the function of SLCGskipA.
 *
 * \tparam n number of generators
 */
template<int n>
struct SLCGskipAt {
 	//! \return modified parameter \c A
	static const unsigned long long a = PowHelper<KMC_SLCG_A, n, n%2==0>::val;
};

template<int n>
struct SLCGskipCtHelper {
	static const unsigned long long val = KMC_SLCG_C*SLCGskipAt<n>::a;
};
/*! \brief Template with the function of SLCGskipC.
 *
 * \tparam n number of generators
 */
template<int n>
struct SLCGskipCt {
 	//! \return modified parameter \c C
	static const unsigned long long c = SumHelper<SLCGskipCtHelper, 0, n, n%2==0, true>::val;
};

/*! \param i internal state of the 64-nit LCG
 * \return next random number in sequence
 */
inline static unsigned long long SLCGen (unsigned long long i) {
	return i*KMC_SLCG_A + KMC_SLCG_C;
}
#include "helper/slcgInst.h"
#else
#pragma KMC_OCL_RUNTIME_MOD paste KMC_SLCG_SKIP
#endif

#endif
