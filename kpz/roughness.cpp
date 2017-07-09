/***************************************************************************
*   Copyright 2011 - 2016 Jeffrey Kelling <j.kelling@hzdr.de>
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

#include "schedulerService.h"

Kpz::Roughness::Roughness(const Roughness* r, int n)
	: w2(0.)
{
	for(auto& a : h)
		a=0.;
	for(int a = 0; a < n; ++a)
	{
		for(int s = 0; s < H_ORDERS; ++s)
			h[s] += r[a].h[s];
		w2 += r[a].w2;
	}
	(*this) /= n;
}

void Kpz::Roughness::normAndSetW2(unsigned long long N)
{
	for(auto& a : h)
		a /= N;
	w2 = h[1] - h[0]*h[0];
}

const Kpz::Roughness& Kpz::Roughness::operator+=(const Kpz::Roughness& other)
{
	for(int i = 0; i < H_ORDERS; ++i)
		h[i] += other.h[i];
	w2 += other.w2;
}

const Kpz::Roughness& Kpz::Roughness::operator/=(double div)
{
	for(auto& a : h)
		a /= div;
	w2 /= div;
}
