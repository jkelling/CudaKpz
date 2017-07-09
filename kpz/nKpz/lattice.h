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

#ifndef N_KPZ_LATTICE_H
#define N_KPZ_LATTICE_H

#include "dim.h"

#include <helper/kmcTemplateMetaHelper.h>
#include <kmcRandom.h>

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>

namespace nKPZ
{

template<int N, int n, class Dim> struct LatticePoint_indexIntern;

template<int N, class Dim>
class LatticePoint
{
	public:

	int c[N]; // coordinates

		LatticePoint() {}
		LatticePoint(int v) {for(int a = 0; a < N; ++a) c[a]=v;}

	inline size_t index(const Dim& dim) const {return LatticePoint_indexIntern<N, N-1, Dim>::fkt(dim, *this);}

		LatticePoint (const Dim dim, size_t index) {
			typename Dim::Mod mod(dim);
			for(int a = N-1; a >= 0; --a) {
				c[a] = mod(index);
				index /= dim;
			}
		}
};

template<int N, class Dim>
std::ostream& operator<<(std::ostream& o, const LatticePoint<N, Dim>& l);

template<int BitsPerSite>
class Pack
{
	public:
	typedef typename TemplateCond<(BitsPerSite < 9), unsigned int, unsigned long long>::Result Data_t;

	protected:
	Data_t m_data;

	public:
	static const int DATA_BITS = (sizeof(Data_t)<<3);
	static const int SITES = DATA_BITS/BitsPerSite;
	static const Data_t M_SITE = (((Data_t)-1)<<(DATA_BITS-BitsPerSite))>>(DATA_BITS-BitsPerSite);
	static const int BITS_PER_SITE = BitsPerSite;

	inline Data_t get(const int site) const {return (m_data>>(site*BITS_PER_SITE))&M_SITE;}
	inline Data_t getBit(const int site, const int n) const {return (m_data>>(site*BITS_PER_SITE+n))&1;}

	/* assume all but the lowest BITS_PER_SITE bits of val to be null */
	inline void set(const int site, const Data_t val) {
		m_data &= ~(M_SITE<<(site*BITS_PER_SITE));
		m_data |= val<<(site*BITS_PER_SITE);
	}
	inline void set(const int site) {
		m_data |= (M_SITE<<(site*BITS_PER_SITE));
	}
	inline int setBit(const int site, const int n) {
		m_data |= 1<<(site*BITS_PER_SITE+n);
	}
	inline void null(const int site) {
		m_data &= ~(M_SITE<<(site*BITS_PER_SITE));
	}
	inline int nullBit(const int site, const int n) {
		m_data &= ~(1<<(site*BITS_PER_SITE+n));
	}

	// next two will blow on overflow
	inline void inc(const int site) {
		const register int t = ((m_data>>(site*BITS_PER_SITE))&M_SITE) + 1;
		m_data &= ~(M_SITE<<(site*BITS_PER_SITE));
		m_data |= t<<(site*BITS_PER_SITE);
	}
	inline void dec(const int site) {
		const register int t = ((m_data>>(site*BITS_PER_SITE))&M_SITE) - 1;
		m_data &= ~(M_SITE<<(site*BITS_PER_SITE));
		m_data |= t<<(site*BITS_PER_SITE);
	}
};

template<int N, class Elem, class Dim>
class Map
{
	protected:
	Elem* m_system;
	size_t m_dim[N];
	Dim m_dimL;
	const int m_upperBound;

	size_t incCoord(const LatticePoint<N, Dim>& l, const int n, const size_t index) {
		if(l.c[n] < m_upperBound)
			return index + m_dim[n];
		else
		{
			LatticePoint<N, Dim> local = l;
			local.c[n] = 0;
			return local.index(m_dimL);
		}
	}
	size_t decCoord(const LatticePoint<N, Dim>& l, const int n, const size_t index) {
		if(l.c[n] > 0)
			return index - m_dim[n];
		else
		{
			LatticePoint<N, Dim> local = l;
			local.c[n] = m_upperBound;
			return local.index(m_dimL);
		}
	}

	public:

		Map(const Dim& dim) : m_system(0), m_dimL(dim), m_upperBound(((int)dim)-1) {
			m_dim[N-1] = 1;
			for(int a = N-2; a >= 0; --a)
				m_dim[a] = dim*m_dim[a+1];
		}
		~Map() {
			delete[] m_system;
		}

	inline void alloc() {m_system = new Elem[size()];}

	inline float& get(const LatticePoint<N, Dim>& l) {return m_system[l.index(m_dimL)];}
	inline float get(const LatticePoint<N, Dim>& l) const {return m_system[l.index(m_dimL)];}
	inline float& get(size_t index) {return m_system[index];}
	inline float get(size_t index) const {return m_system[index];}
	
	inline size_t dim(int n) const {return m_dim[n];}
	inline Dim dim() const {return m_dimL;}
	inline size_t size() const {return m_dimL*m_dim[0];}

	inline void randomizePoint(LatticePoint<N, Dim>& l) {
		for(int a = 0; a < N; ++a)
		{
			l.c[a] = (dsfmt_gv_genrand_close_open() * m_dimL);
		}
	}

};

template<int N, int BitsPerSite, class Dim>
class PackMap : public Map <N, Pack<BitsPerSite>, Dim>
{
	public:
	static const int BITS_PER_SITE = BitsPerSite;

	protected:
	using Map<N, Pack<BITS_PER_SITE>, Dim>::m_system;
	using Map<N, Pack<BITS_PER_SITE>, Dim>::m_dimL;

	public:

	using Map<N, Pack<BITS_PER_SITE>, Dim>::size;

		PackMap(const Dim& dim) : Map<N, Pack<BITS_PER_SITE>, Dim>(dim) {
			if(!dim.isEven()) // dim has to be even
				throw std::runtime_error("FATAL: Lateral system dimension has to be even.");
			size_t s = size()/Pack<BITS_PER_SITE>::SITES;
			if(size() % Pack<BITS_PER_SITE>::SITES)
				++s;
			m_system = new Pack<BITS_PER_SITE>[s];
			memset(m_system, 0, s*sizeof(Pack<BITS_PER_SITE>));
		}

	inline int get(const size_t index) const {
		return m_system[index/Pack<BITS_PER_SITE>::SITES].get(index%Pack<BITS_PER_SITE>::SITES);
	}
	inline int get(const LatticePoint<N, Dim>& l) const {return get(l.index(m_dimL));}
	inline void set(const size_t index, const int val) {
		m_system[index/Pack<BITS_PER_SITE>::SITES].set(index%Pack<BITS_PER_SITE>::SITES, val);
	}
	inline void set(const size_t index) {
		m_system[index/Pack<BITS_PER_SITE>::SITES].null(index%Pack<BITS_PER_SITE>::SITES);
	}
	inline void null(const size_t index) {
		m_system[index/Pack<BITS_PER_SITE>::SITES].null(index%Pack<BITS_PER_SITE>::SITES);
	}

	std::ostream& printMemStat(std::ostream& o);
};

#include "lattice.cpp"
}

#endif
