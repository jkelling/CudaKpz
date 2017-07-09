/***************************************************************************
*   Copyright 2014 Jeffrey Kelling <j.kelling@hzdr.de>
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

#include "CPURngIface.h"

#include <tinyMT/tinymt32_param.h>

#include <cstdlib>
#include <climits>

Kpz::RngIface::DSFMT* Kpz::RngIface::DSFMT::create(int n, dsfmt_t* dsfmt)
{
	DSFMT* ret = new DSFMT[n];
	for(int a = 0; a < n; ++a)
		ret[a].seed(dsfmt_genrand_close_open(dsfmt)*INT_MAX);
	return ret;
}

void Kpz::RngIface::TinyMT32::setMat(int index)
{
	m_tmt.mat1 = tinyMTmat1[index];
	m_tmt.mat2 = tinyMTmat2[index];
	m_tmt.tmat = tinyMTtmat[index];
}
