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

#ifndef N_KPZ_SLOPE_SIM_H
#define N_KPZ_SLOPE_SIM_H

#include "lattice.h"

namespace nKPZ
{

template<int N, class Dim>
class SlopeSystem : public PackMap<N, N, Dim>
{
	public:
	using PackMap<N, N, Dim>::BITS_PER_SITE;

	private:
	template<int n>
	inline void incNeighbours(LatticePoint<N, Dim>& l, const size_t index);

	using Map<N, Pack<BITS_PER_SITE>, Dim>::m_system;
	using Map<N, Pack<BITS_PER_SITE>, Dim>::m_dim;
	using Map<N, Pack<BITS_PER_SITE>, Dim>::m_dimL;
	using Map<N, Pack<BITS_PER_SITE>, Dim>::m_upperBound;

	public:

	using Map<N, Pack<BITS_PER_SITE>, Dim>::randomizePoint;
	using Map<N, Pack<BITS_PER_SITE>, Dim>::size;
	using PackMap<N, BITS_PER_SITE, Dim>::get;
	using PackMap<N, BITS_PER_SITE, Dim>::set;
	using PackMap<N, BITS_PER_SITE, Dim>::null;

	inline int getBit(const size_t index, const int n) const {
		return m_system[index/Pack<BITS_PER_SITE>::SITES].getBit(index%Pack<BITS_PER_SITE>::SITES, n);
	}
	inline int getBit(const LatticePoint<N, Dim>& l, const int n) const {return getBit(l.index(m_dimL), n);}
	inline int setBit(const size_t index, const int n) {
		return m_system[index/Pack<BITS_PER_SITE>::SITES].setBit(index%Pack<BITS_PER_SITE>::SITES, n);
	}
	inline int nullBit(const size_t index, const int n) {
		return m_system[index/Pack<BITS_PER_SITE>::SITES].nullBit(index%Pack<BITS_PER_SITE>::SITES, n);
	}

		SlopeSystem(const Dim& dim) : PackMap<N, BITS_PER_SITE, Dim>(dim) {}

	void update();
	void init(); 
	double roughness() const;

	template<int, int, class> friend struct SlopeSystem_setNeighboursIterate;
};

#include "slopeSim.cpp"
}
#endif
