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

#include "KMCh5.h"

#include "dSFMT/dSFMT.h"

#include <H5Cpp.h>

#include <iostream>
#include <memory>

static std::string mkH5DsName(const std::string& prefix)
{
	if(!prefix.empty())
		return prefix + ".dsfmt";
	else
		return "dsfmt";
}

void Kmc::writeDSFMT(H5::CommonFG& dest, dsfmt_t* dsfmt, const std::string& prefix)
{
	const auto name = mkH5DsName(prefix);

	const hsize_t size = sizeof(dsfmt_t)/4;
	H5::DataSpace dspace(1, &size);
	auto dset = dest.createDataSet(name, H5::PredType::STD_U32LE, dspace);
	dset.write(dsfmt, H5::PredType::STD_U32LE);
}

void Kmc::readDSFMT(H5::CommonFG& src, dsfmt_t* dsfmt, const std::string& prefix)
{
	const auto name = mkH5DsName(prefix);

	auto dset = src.openDataSet(name);
	auto dspace = dset.getSpace();
	if(dspace.getSimpleExtentNdims() != 1 )
		throw std::domain_error("[readDSFMT] Invalid dsfmt DataSet: Wrong NDim");
	hsize_t size;
	dspace.getSimpleExtentDims(&size);
	if(size != sizeof(dsfmt_t)/4)
		throw std::range_error("[readDSFMT] Invalid dsfmt state in HDF5 file (wrong size).\n");

	dset.read(dsfmt, H5::PredType::STD_U32LE);

	std::cout << "[readDSFMT] restored DSFMT state from file.\n";
}
