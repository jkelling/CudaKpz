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

template<int N, class Dim> template<int n>
inline void System<N, Dim>::incNeighbours(LatticePoint<N, Dim>& l, const size_t index) {
	if(l.c[n] < m_upperBound)
		inc(index + m_dim[n]);
	else
	{
		l.c[n] = 0;
		inc(l.index(m_dim[N-2]));
		l.c[n] = m_upperBound;
	}
	if(l.c[n] > 0)
		inc(index - m_dim[n]);
	else
	{
		l.c[n] = m_upperBound;
		inc(l.index(m_dim[N-2]));
		l.c[n] = 0;
	}
}
template<int N, int n, class Dim> struct System_incNeighboursIterate {
	static inline void fkt(LatticePoint<N, Dim>& l, const size_t index, System<N, Dim>& s) {
		s.template incNeighbours<n>(l, index);
		System_incNeighboursIterate<N, n-1, Dim>::fkt(l, index, s);
	}
};
template<int N, class Dim> struct System_incNeighboursIterate<N, 0, Dim> {
	static inline void fkt(LatticePoint<N, Dim>& l, const size_t index, System<N, Dim>& s) {
		s.template incNeighbours<0>(l, index);
	}
};

template<int N, int n, class Dim> struct System_init {
	static inline void fkt(LatticePoint<N, Dim> l, int sign, System<N, Dim>& s) {
		for(; l.c[n] < (int)s.dim(); ++l.c[n], ++sign)
			System_init<N, n-1, Dim>::fkt(l, sign&1, s);
	}
};
template<int N, class Dim> struct System_init<N,0,Dim> {
	static inline void fkt(LatticePoint<N,Dim> l, int sign, System<N,Dim>& s) {
		static const int val[] = {0,2*N};
		for(; l.c[0] < (int)s.dim(); ++l.c[0], ++sign)
			s.set(l.index(s.dim()), val[sign&1]);
	}
};
