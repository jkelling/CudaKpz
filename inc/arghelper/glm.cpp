/***************************************************************************
*   Copyright (C) 2013 Jeffrey Kelling <kelling.jeffrey@ages-skripte.org>
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

#include "glm.h"

int getArg(glm::dvec3& val, int i, int argc, char* argv[], int flags)
{
	for(int a = 0; a < 3; ++a)
	{
		if(!getArg(val[a], i, argc, argv, flags))
			return a;
		++i;
	}
	return 3;
}
