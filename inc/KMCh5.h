/***************************************************************************
*   Copyright 2016 Jeffrey Kelling <j.kelling@hzdr.de>
*                  Helmholtz-Zentrum Dresden-Rossendorf
*                  Department of Information Services and Computing
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

#pragma once

#include <string>

struct DSFMT_T;
typedef DSFMT_T dsfmt_t;

namespace H5 {
	class CommonFG;
}

namespace Kmc
{
	void writeDSFMT(H5::CommonFG& dest, dsfmt_t* dsfmt, const std::string& prefix = "");
	void readDSFMT(H5::CommonFG& src, dsfmt_t* dsfmt, const std::string& prefix = "");
}
