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

template<int N, int n, class Dim> struct SlopeSystem_setNeighboursIterate {
	static inline void fkt(LatticePoint<N, Dim>& l, const size_t index, SlopeSystem<N, Dim>& s) {
		s.setBit(s.decCoord(l, n, index),n);
		SlopeSystem_setNeighboursIterate<N, n-1, Dim>::fkt(l, index, s);
	}
};
template<int N, class Dim> struct SlopeSystem_setNeighboursIterate<N, 0, Dim> {
	static inline void fkt(LatticePoint<N, Dim>& l, const size_t index, SlopeSystem<N, Dim>& s) {
		s.setBit(s.decCoord(l, 0, index),0);
	}
};

template<int N, class Dim>
void SlopeSystem<N, Dim>::update() {
	LatticePoint<N, Dim> l;
	randomizePoint(l);
	size_t index = l.index(m_dimL);
	
	if(get(index) != Pack<BITS_PER_SITE>::M_SITE)
		return;

	for(int a = N-1; a >= 0; --a)
		if(getBit(decCoord(l, a, index), a))
			return;

	null(index);
	SlopeSystem_setNeighboursIterate<N, N-1, Dim>::fkt(l, index, *this);
}

template<int N, int n, class Dim> struct SlopeSystem_init {
	static inline void fkt(LatticePoint<N, Dim> l, int sign, SlopeSystem<N, Dim>& s) {
		for(; l.c[n] < (int)s.dim(); ++l.c[n], ++sign)
			SlopeSystem_init<N, n-1, Dim>::fkt(l, sign&1, s);
	}
};
template<int N, class Dim> struct SlopeSystem_init<N,0,Dim> {
	static inline void fkt(LatticePoint<N,Dim> l, int sign, SlopeSystem<N,Dim>& s) {
		static const int val[] = {0,Pack<SlopeSystem<N,Dim>::BITS_PER_SITE>::M_SITE};
		for(; l.c[0] < (int)s.dim(); ++l.c[0], ++sign)
			s.set(l.index(s.dim()), val[sign&1]);
	}
};
template<int N, class Dim>
void SlopeSystem<N, Dim>::init()
{
	SlopeSystem_init<N, N-1, Dim>::fkt(LatticePoint<N, Dim>(0), 0, *this);
}

template<int N, int n, class Dim> struct SlopeSystem_roughness {
	static inline void fkt(LatticePoint<N, Dim> l, double& sum, double& sqsum, const SlopeSystem<N, Dim>& s, int h = 0) {
		for(; l.c[n] < (int)s.dim(); ++l.c[n])
		{
			//std::cerr << l << ' ' << h << ' ' << s.getBit(l, n) <<  '\n';
			SlopeSystem_roughness<N, n-1, Dim>::fkt(l, sum, sqsum, s, h);
			h += (s.getBit(l, n) ? +1 : -1);
		}
	}
};
template<int N, class Dim> struct SlopeSystem_roughness<N,0,Dim> {
	static inline void fkt(LatticePoint<N, Dim> l, double& sum, double& sqsum, const SlopeSystem<N, Dim>& s, int h = 0) {
		static const int val[] = {0,Pack<SlopeSystem<N, Dim>::BITS_PER_SITE>::M_SITE};
		for(; l.c[0] < (int)s.dim(); ++l.c[0])
		{
			const double v = h;
			sum += v;
			sqsum += v*v;
			h += (s.getBit(l, 0) ? +1 : -1);
		}
	}
};
template<int N, class Dim>
double SlopeSystem<N, Dim>::roughness() const
{
	double sum = 0, sqsum = 0;
	SlopeSystem_roughness<N, N-1, Dim>::fkt(LatticePoint<N, Dim>(0), sum, sqsum, *this);
	sum /= size();
	sqsum /= size();
	return (sqsum-sum*sum);
}
