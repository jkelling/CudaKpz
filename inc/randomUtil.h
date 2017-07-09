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

/*! \file randomUtil.h
 * \brief Specifically distributed random numbers.
 * \ingroup DGutils
 */

#ifndef KMC_RANDOM_UTIL_H
#define KMC_RANDOM_UTIL_H

#include <dSFMT/dSFMT.h>

#include <cmath>

#ifndef EXP_RAND_MAX_TRIES
#define EXP_RAND_MAX_TRIES 1
#endif
const double EXP_RAND_MAX = 10.;
/*! \tparam Rand Random number generator predicate.
 * \param lam decay parameter
 * \param rand instance of \p Rand
 * \param max maximum value for the distribution
 * \returns random number distributed accord	ing to \f$\exp(-x/\lambda)\f$
 */
template<class Rand>
double expRand(double lam, Rand rand, double max = EXP_RAND_MAX)
{
	double r;
	for(int a = EXP_RAND_MAX_TRIES; a > 0; --a)
	{
		r = rand();
		r = -log(r)*lam;
		if(r < max)
			return r;
	}
	// screw it, this small error should be ok, adjust EXP_RAND_MAX_TRIES above if not
	// is there an error at all? since exp is self congruent...
	return fmod(r, max); 
}
/*! expRand() using dSFMT.
 * \param lam decay parameter
 * \param max maximum value for the distribution
 * \returns random number distributed according to \f$\exp(-x/\lambda)\f$
 */
inline double expRandDSFMT(double lam, double max = EXP_RAND_MAX)
{
	return expRand(lam, dsfmt_gv_genrand_open_close, max);
}

/*! \tparam Rand Random number generator predicate.
 * \param lam decay parameter
 * \param rand instance of \p Rand
 * \param max maximum value for the distribution
 * \returns random number distributed according to \f$\exp(-x/\lambda)\f$, where
 * values larger than \p max will be mapped to `max-0`
 */
template<class Rand>
double expRandHardSat(double lam, Rand rand, double max = EXP_RAND_MAX)
{
	double r = rand();
	r = -log(r)*lam;
	if(r >= max)
		r = max*.999999;
	return r;
}
/*! expRandHardSat() using dSFMT.
 * \param lam decay parameter
 * \param max maximum value for the distribution
 * \returns random number distributed according to \f$\exp(-x/\lambda)\f$
 */
inline double expRandHardSatDSFMT(double lam, double max = EXP_RAND_MAX)
{
	return expRandHardSat(lam, dsfmt_gv_genrand_open_close, max);
}
#endif
