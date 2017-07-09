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

#include "KMCsplash.h"

#include "dSFMT/dSFMT.h"

#include <iostream>

// ::splash::ColTypeUInt32 Kmc::splash::ColTypeUInt32;

void KmcSplash::writeDSFMT(splash::DataCollector* data, int id, dsfmt_t* dsfmt, const std::string& prefix)
{
	std::string name = prefix + ".dsfmt";
	data->write(id, ColTypeUInt32, 1, splash::Selection(splash::Dimensions(sizeof(dsfmt_t)/4,1,1)), name.c_str(), dsfmt);
}

bool KmcSplash::readDSFMT(splash::DataCollector* data, int id, dsfmt_t* dsfmt, const std::string& prefix)
{
	std::string name = prefix + ".dsfmt";
	splash::Dimensions size;
	try {
		data->read(id, name.c_str(), size, 0);
	}
	catch (splash::DCException) {
		std::cerr << "[readDSFMT_H5][EE] Data does not provide dsfmt at id " << id << "\n";
		return false;
	}
	if(size[0]*4 != sizeof(dsfmt_t) || size[1] > 1 || size[2] > 1)
	{
		std::cerr << "[readDSFMT_H5][EE] stored dsfmt has incompatible size.\n";
		return false;
	}
	data->read(id, name.c_str(), size, dsfmt);
	return true;
}

std::string KmcSplash::getFullSerialFilename(unsigned int rankX, const std::string& baseFilename)
{
	std::stringstream serial_filename;
	serial_filename << baseFilename << "_" << rankX << "_0_0.h5";

	return serial_filename.str();
}
