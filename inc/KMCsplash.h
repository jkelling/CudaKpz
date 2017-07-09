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

#ifndef KMC_SPLASH_H
#define KMC_SPLASH_H

#include <splash/splash.h>
#include <string>

struct DSFMT_T;
typedef DSFMT_T dsfmt_t;

namespace KmcSplash
{
	const splash::ColTypeUInt32 ColTypeUInt32;
	const splash::ColTypeDouble ColTypeDouble;

	void writeDSFMT(splash::DataCollector* data, int id, dsfmt_t* dsfmt, const std::string& prefix = "");
	bool readDSFMT(splash::DataCollector* data, int id, dsfmt_t* dsfmt, const std::string& prefix = "");

	std::string getFullSerialFilename(unsigned int rankX, const std::string& baseFilename);
}
#endif
