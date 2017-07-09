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

template<int N, int n, class Dim> struct LatticePoint_indexIntern {
	static inline size_t fkt (const Dim dim, const LatticePoint<N, Dim>& l) {
		return LatticePoint_indexIntern<N, n-1, Dim>::fkt(dim, l)*dim + l.c[n];
	}
};
template<int N, class Dim> struct LatticePoint_indexIntern<N, 0, Dim> {
	static inline size_t fkt(const Dim dim, const LatticePoint<N, Dim>& l) {
		return l.c[0];
	}
};

template<int N, class Dim>
std::ostream& operator<<(std::ostream& o, const LatticePoint<N, Dim>& l)
{
	int a = N-1;
	o << "( " << l.c[a];
	for(--a; a >= 0; --a)
		o << ' ' << l.c[a];
	return o << " )";
}

template<int N, int BitsPerSite, class Dim>
std::ostream& PackMap<N, BitsPerSite, Dim>::printMemStat(std::ostream& o)
{
	return o << "System has " << size() << " sites (" << BITS_PER_SITE << " bits per site), using " 
		<< size()/Pack<BITS_PER_SITE>::SITES*sizeof(typename Pack<BITS_PER_SITE>::Data_t) << " bytes.\n";
}
