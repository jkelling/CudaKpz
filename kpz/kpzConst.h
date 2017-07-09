/***************************************************************************
*   Copyright 2011 - 2014 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef KPZ_CONST_H
#define KPZ_CONST_H

// single hit double tiling at block level
#define KPZ_INNER_DT 100
// the following are only implemented fo GPU DT
// single hit double tiling (balanced) dead border at block level
#define KPZ_INNER_DTDBSH 101
// dead border (fixed) with in-sequence delayed borders at block level
#define KPZ_INNER_DB 200
// single-hit dead border with in-sequence delayed at block level
#define KPZ_INNER_DBSH 201
// single-hit dead border with out-of-sequence replaced borders at block level
#define KPZ_INNER_DBSHREP 211

#include "param.h"

#define KMC_NO_RNG_SHUFFLING
#define KMC_L_LCG_N 4
#include <kmcRandom.h>

#ifndef __OPENCL_VERSION__
namespace Kpz
{
	//all dimensions are those of the slope-grid
	#ifdef KPZ_FERMI
		const int L_SCA_MAX_THREADS = 11;
	#elif defined KPZ_K80
		const int L_SCA_MAX_THREADS = 11;
	#else
		const int L_SCA_MAX_THREADS = 9;
	#endif
	const int L_BASIC_CELL_DIM_X = 2; // double word, 2 bits per site
	const int BASIC_CELL_DIM_X = 1<<L_BASIC_CELL_DIM_X;
	const int L_BASIC_CELL_DIM_Y = 2;
	const int BASIC_CELL_DIM_Y = 1<<L_BASIC_CELL_DIM_Y;

	const int DEAD_BORDER_DIM = 1;

#if KPZ_INNER_DD == KPZ_INNER_DT
	static constexpr const char* INNER_DD_STRING = "DT";
	static constexpr const char* INNER_DD_STRING_LONG = "single-hit double tiling";
#elif KPZ_INNER_DD == KPZ_INNER_DTDBSH
	static constexpr const char* INNER_DD_STRING = "DTDBSH";
	static constexpr const char* INNER_DD_STRING_LONG = "single-hit double tiling (balanced) dead border";
#elif KPZ_INNER_DD == KPZ_INNER_DB
	static constexpr const char* INNER_DD_STRING = "DB/NB";
	static constexpr const char* INNER_DD_STRING_LONG = "delayed border / no border if p=1,q=0";
#elif KPZ_INNER_DD == KPZ_INNER_DBSH
	static constexpr const char* INNER_DD_STRING = "DBSH/NB";
	static constexpr const char* INNER_DD_STRING_LONG = "single-hit delayed border / no border if p=1,q=0";
#elif KPZ_INNER_DD == KPZ_INNER_DBSHREP
	static constexpr const char* INNER_DD_STRING = "DBSHRep/NB";
	static constexpr const char* INNER_DD_STRING_LONG = "single-hit delayed border, borders out-of-sequence / no border if p=1,q=0";
#else
#error "KPZ_INNER_DD is invalid."
#endif
}
#endif

typedef unsigned long long KmcRandom_t;
const int KMC_MAX_BLOCKS = 30; // max number of blocks expected (Tesla), see KMC_BOCKS
// const int KPZ_SLCG_SKIP = KMC_MAX_BLOCKS<<Kpz::L_THREADS;
const int KPZ_SLCG_SKIP = 30721; //! good for up to 30 blocks each 1024 threads
const KmcRandom_t KMC_SLCG_A_SKIP = SLCGskipAt<KPZ_SLCG_SKIP>::a;
const KmcRandom_t KMC_SLCG_C_SKIP = SLCGskipCt<KPZ_SLCG_SKIP>::c;

#endif
