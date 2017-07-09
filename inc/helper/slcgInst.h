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


// Instantiation for KmcDBConst<3,2,2>
template<>
struct SLCGskipAt<3840> {
	static const unsigned long long a = 12972124132769090561ull;
};
template<>
struct SLCGskipCt<3840> {
	static const unsigned long long c = 1151287398592190208ull;
};

// Instantiation for kmcConst 3 3 2
template<>
struct SLCGskipAt<7680> {
	static const unsigned long long a = 8477683129405052929ull;
};
template<>
struct SLCGskipCt<7680> {
	static const unsigned long long c = 11710633003335364096ull;
};

// Instatiation for kmcConst 3 3 3 (Fermi), kpz 5 4
template<>
struct SLCGskipAt<15360> {
	static const unsigned long long a = 4707689352879484929ull;
};

template<>
struct SLCGskipCt<15360> {
	static const unsigned long long c = 4920721155844365312ull;
};

// Instatiation for kpz 5 5 (Fermi)
template<>
struct SLCGskipAt<30720> {
	static const unsigned long long a = 4602270209800822785ull;
};

template<>
struct SLCGskipCt<30720> {
	static const unsigned long long c = 8580373701247416320ull;
};

// Instatiation for kpz 6 5 (GK210)
template<>
struct SLCGskipAt<61440> {
	static const unsigned long long a = 11348253892226105345ull;
};

template<>
struct SLCGskipCt<61440> {
	static const unsigned long long c = 1524031857202130944ull;
};

// Instatiation for kpz 5 5 (Fermi) ... odd
template<>
struct SLCGskipAt<30721> {
	static const unsigned long long a = 16320982315376414973ull;
};

template<>
struct SLCGskipCt<30721> {
	static const unsigned long long c = 11667173346849256211ull;
};
