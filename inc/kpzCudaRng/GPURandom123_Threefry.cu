/***************************************************************************
*   Copyright 2014 - 2016 Jeffrey Kelling <j.kelling@hzdr.de>
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

#include "GPURandom123_Threefry.h"
#include "GPURandom123_Threefry.cuh"

#include <cassert>
#include <climits>
#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <stdexcept>
#include <vector>

#include <cuda.h>
#include <cudaError.h>

/*#include <H5Cpp.h>*/

void Random123_Threefry::initCUDA()
{
	{ // round up to multiples or warps
		int tmp = 0;
		cudaDeviceProp prop;
		CUDA_SAFE_CALL_THROW(cudaGetDevice(&tmp));
		CUDA_SAFE_CALL_THROW(cudaGetDeviceProperties(&prop, tmp));
		if(m_generatorCount%prop.warpSize)
			m_generatorCount = (m_generatorCount/prop.warpSize+1)*prop.warpSize;
	}

	CUDA_SAFE_CALL(cudaMalloc((void**)&d_Random, m_generatorCount*sizeof(unsigned int)));
	if(!d_Random)
	{
		throw std::runtime_error("Insufficient devicememory for Random123_Threefry state.");
	}
	CUDA_SAFE_CALL_THROW( cudaMemset(d_Random, 0, m_generatorCount*sizeof(unsigned int)) );
}

void Random123_Threefry::randomize()
{
	assert(m_dsfmt != 0);
	m_key[0] = dsfmt_genrand_close_open(m_dsfmt)*UINT_MAX;
	m_key[1] = dsfmt_genrand_close_open(m_dsfmt)*UINT_MAX;
}

bool Random123_Threefry::minGenerators(unsigned int ngenmin)
{
	if(m_generatorCount >= ngenmin)
		return true;

	if (d_Random)
		CUDA_SAFE_CALL(cudaFree(d_Random));
	m_generatorCount = ngenmin;
	initCUDA();
	return false;
}

Random123_Threefry::~Random123_Threefry()
{
	if (d_Random)
		CUDA_SAFE_CALL(cudaFree(d_Random));
}

#if 0
#ifdef USE_LIB_SPLASH
#include <KMCsplash.h>

void Random123_Threefry::writeH5(splash::DataCollector* data, int id, const char* prefix)
{
	std::string name;
	if(prefix)
		name = prefix;
	name += ".tinyMT";

	auto rnd = std::unique_ptr<unsigned int[]>( new unsigned int[m_generatorCount*7] );
	CUDA_SAFE_CALL_THROW( cudaMemcpy(rnd.get(), d_Random, m_generatorCount*sizeof(unsigned int)*7, cudaMemcpyDeviceToHost) );
	data->write(id, KmcSplash::ColTypeUInt32, 2, splash::Selection(splash::Dimensions(m_generatorCount, 7, 1)), name.c_str(), rnd.get());
}

bool Random123_Threefry::readH5(splash::DataCollector* data, int id, const char* prefix)
{
	std::string name;
	if(prefix)
		name = prefix;
	name += ".tinyMT";

	splash::Dimensions size;
	try {
		data->read(id, name.c_str(), size, 0);
	}
	catch (splash::DCException) {
		return false;
	}

	if(size[0] != m_generatorCount)
	{
		std::cerr << "[GPUTinyMT][readH5] incompatible tinyMT state in HDF5 file.\n";
		return false;
	}
	auto rnd = std::unique_ptr<unsigned int[]>( new unsigned int[m_generatorCount*7] );
	data->read(id, name.c_str(), size, rnd.get());
	CUDA_SAFE_CALL_THROW( cudaMemcpy(d_Random, rnd.get(), m_generatorCount*sizeof(unsigned int)*7, cudaMemcpyHostToDevice) );
	std::cout << "[GPUTinyMT][readH5] restored TinyMT state from file.\n";
	return true;
}
#endif

#include <H5Cpp.h>

static std::string mkH5DsName(const char* prefix)
{
	std::string name;
	if(prefix)
	{
		name = prefix;
		name += ".tinyMT";
	}
	else
		name = "tinyMT";

	return name;
}

void TinyMT32::writeH5(H5::CommonFG& dest, const char* prefix)
{
	const auto name = ::mkH5DsName(prefix);

	const hsize_t size[2] = {m_generatorCount, 7};
	auto rnd = std::unique_ptr<unsigned int[]>( new unsigned int[size[0]*size[1]] );
	CUDA_SAFE_CALL_THROW( cudaMemcpy(rnd.get(), d_Random, size[0]*size[1]*sizeof(unsigned int), cudaMemcpyDeviceToHost) );

	H5::DataSpace dspace(2, size);
	auto dset = dest.createDataSet(name, H5::PredType::STD_U32LE, dspace);
	dset.write(rnd.get(), H5::PredType::STD_U32LE);
}

void TinyMT32::readH5(H5::CommonFG& src, const char* prefix)
{
	const auto name = ::mkH5DsName(prefix);

	auto dset = src.openDataSet(name);
	auto dspace = dset.getSpace();
	if(dspace.getSimpleExtentNdims() != 2 )
		throw std::domain_error("[TinyMT32::readH5] Invalid tinyMT DataSet: Wrong NDim");
	hsize_t size[2];
	dspace.getSimpleExtentDims(size);
	if(size[1] != 7)
		throw std::domain_error("[TinyMT32::readH5] Invalid tinyMT DataSet: Wrong number of sections ([1] != 7).");
	if(size[0] != m_generatorCount)
		throw std::range_error("[TinyMT32::readH5] Incompatible tinyMT state in HDF5 file.\n");

	auto rnd = std::unique_ptr<unsigned int[]>( new unsigned int[m_generatorCount*7] );
	dset.read(rnd.get(), H5::PredType::STD_U32LE);
	CUDA_SAFE_CALL_THROW( cudaMemcpy(d_Random, rnd.get(), m_generatorCount*sizeof(unsigned int)*7, cudaMemcpyHostToDevice) );

	std::cout << "[GPUTinyMT][readH5] restored TinyMT state from file.\n";
}
#endif
