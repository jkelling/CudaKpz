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

#ifndef N_KPZ_SUM_SIM_H
#define N_KPZ_SUM_SIM_H

#include "lattice.h"

namespace nKPZ
{

template<int N, int n, class Dim> struct System_incNeighboursIterate;
template<int N, int n, class Dim> struct System_init;

template<int N, class Dim>
class System : public PackMap<N, Log2RoundUp<N*2>::val, Dim>
{
	public:
	using PackMap<N, Log2RoundUp<N*2>::val, Dim>::BITS_PER_SITE;

	private:
	template<int n>
	inline void incNeighbours(LatticePoint<N, Dim>& l, const size_t index);

	using Map<N, Pack<BITS_PER_SITE>, Dim>::m_system;
	using Map<N, Pack<BITS_PER_SITE>, Dim>::m_dim;
	using Map<N, Pack<BITS_PER_SITE>, Dim>::m_dimL;
	using Map<N, Pack<BITS_PER_SITE>, Dim>::m_upperBound;

	public:

	using Map<N, Pack<BITS_PER_SITE>, Dim>::randomizePoint;
	using PackMap<N, BITS_PER_SITE, Dim>::get;
	using PackMap<N, BITS_PER_SITE, Dim>::set;
	using PackMap<N, BITS_PER_SITE, Dim>::null;

	inline void inc(const size_t index) {
		m_system[index/Pack<BITS_PER_SITE>::SITES].inc(index%Pack<BITS_PER_SITE>::SITES);
	}
	inline void dec(const size_t index) {
		m_system[index/Pack<BITS_PER_SITE>::SITES].dec(index%Pack<BITS_PER_SITE>::SITES);
	}

		System(const Dim& dim) : PackMap<N, BITS_PER_SITE, Dim>(dim) {}

	void update() {
		LatticePoint<N, Dim> l;
		randomizePoint(l);
		size_t index = l.index(m_dimL);
		
		if(get(index) == (N*2))
		{
			null(index);
			System_incNeighboursIterate<N, N-1, Dim>::fkt(l, index, *this);
		}
	}

	void init() {
		System_init<N, N-1, Dim>::fkt(LatticePoint<N, Dim>(0), 0, *this);
	}

	template<int N2, int n, class Dim2> friend struct System_incNeighboursIterate;
};

template<int N, class Dim>
class HeightMap : public Map<N, float, Dim>
{
	using Map<N, float, Dim>::m_system;
	using Map<N, float, Dim>::m_dim;

	void tryNeighbours(const size_t index, const System<N, Dim>& system, const bool min);
	void tryNeighbours(const size_t index, const System<N, Dim>& system);
	void trySet(const size_t index, const System<N, Dim>& system);

	public:

	using Map<N, float, Dim>::dim;
	using Map<N, float, Dim>::size;

		HeightMap(const Dim& dim) : Map<N, float, Dim>(dim) {
			Map<N, float, Dim>::alloc();
			set();
		}

	void set(const float val = NAN) {
		for(size_t a = 0; a < size(); ++a)
			m_system[a] = val;
	}

	void fromSystem(const System<N, Dim>& system) {
		size_t index;
		for(index = 0; index < size(); ++index) {
			const int t = system.get(index);
			if(t == 0 || t == (N*2)) {
				m_system[index] = 0; // arbitrary reference height
				tryNeighbours(index, system, t == 0);
				break;
			}
		}
	}

	long double roughness() {
		long double sum = 0, sumsq = 0;
		for(size_t a = 0; a < size(); ++a) {
			const long double t = m_system[a];
			sum += t;
			sumsq += t*t;
		}
		sumsq /= size();
		sum /= size();
		return (sumsq - sum*sum);
	}

};

#include "sumSim.cpp"
#include "heightMap.cpp"
}
#endif
