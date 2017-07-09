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

#ifndef KPZ_KPZ_SIMULATION_PRIMITIVES
#define KPZ_KPZ_SIMULATION_PRIMITIVES

#include "kpzConst.h"

namespace Kpz
{

// bitselctionmask for koordinatepicking
const int O_MASK = 0xFF;
const float O_RND_SUP = O_MASK+1;
const int O_WIDTH = 8;

__device__ static int shift(int x, int y)
{
	return ((x&3)<<1) // select row by x
		| ((y&3)<<3); // select column by y
}
}
#endif
