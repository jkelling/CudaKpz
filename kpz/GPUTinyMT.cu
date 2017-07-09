/***************************************************************************
*   Copyright 2014 - 2014 Jeffrey Kelling <j.kelling@hzdr.de>
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

#include "GPUTinyMT.h"

#include "kpzRandom.cu"
#include "kpzConst.h"

#include <tinyMT/tinymt32_param.h>

#include <cassert>
#include <climits>
#include <iostream>
#include <string>
#include <memory>

#include <cuda.h>
#include <cudaError.h>

void TinyMT32::initCUDA()
{
	CUDA_SAFE_CALL(cudaMalloc((void**)&d_Random, m_generatorCount*sizeof(unsigned int)*7));
	CUDA_SAFE_CALL(cudaMalloc((void**)&d_tmat, m_generatorCount*sizeof(unsigned int)));
	if(!d_Random || !d_tmat)
	{
		std::cerr << "Insufficient devicememory for TinyMT state.\n";
		exit(1);
	}
	CUDA_SAFE_CALL_THROW( cudaMemcpy(d_Random+4*m_generatorCount, \
				tinyMTmat1, m_generatorCount*sizeof(unsigned int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL_THROW( cudaMemcpy(d_Random+5*m_generatorCount, \
				tinyMTmat2, m_generatorCount*sizeof(unsigned int), cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL_THROW( cudaMemcpy(d_Random+6*m_generatorCount, \
				tinyMTtmat, m_generatorCount*sizeof(unsigned int), cudaMemcpyHostToDevice) );
	std::cout << "TinyMT32: Using parametersets: 0 .. " << m_generatorCount-1 << '\n';

	int device;
	CUDA_SAFE_CALL(cudaGetDevice(&device));
	cudaDeviceProp prop;
	CUDA_SAFE_CALL(cudaGetDeviceProperties(&prop, device));
	m_blocks = prop.multiProcessorCount;
	m_threads = m_generatorCount/m_blocks;
	if(m_generatorCount%m_blocks)
		m_blocks += 1;
}

__global__ void seedGPU(unsigned int* d_Random, size_t generatorCount)
{
	if(blockIdx.x*blockDim.x+threadIdx.x < generatorCount)
	{
		TinyMT32::Device rng(d_Random);
			rng.seed();
	}
}

void TinyMT32::randomize()
{
	assert(m_generatorCount <= KMC_N_TMT_PARAM_SETS);
	auto rnd = std::unique_ptr<unsigned int[]>( new unsigned int[m_generatorCount] );
	for (int a = 0; a<m_generatorCount; ++a)
	{
		rnd[a] = dsfmt_genrand_close_open(m_dsfmt)*UINT_MAX;
	}
	CUDA_SAFE_CALL_THROW( cudaMemcpy(d_Random, rnd.get(), m_generatorCount*sizeof(unsigned int), cudaMemcpyHostToDevice) );

	// determine optimal <<<blocks,threads>>> for seeding
	cudaDeviceProp prop;
	int tmp;
	CUDA_SAFE_CALL_THROW(cudaGetDevice(&tmp));
	CUDA_SAFE_CALL_THROW(cudaGetDeviceProperties(&prop, tmp));
	int threads = prop.maxThreadsPerBlock;
	int blocks = m_generatorCount/threads;
	if(m_generatorCount%threads)
		++blocks;
	if(blocks < prop.multiProcessorCount)
		blocks = prop.multiProcessorCount;
	threads = m_generatorCount/blocks;
	if(m_generatorCount%blocks)
		++threads;
	if(threads%prop.warpSize)
		threads = (threads/prop.warpSize+1)*prop.warpSize;

	seedGPU<<<blocks, threads>>> (d_Random, m_generatorCount);
	CUDA_SAFE_CALL_THROW( cudaGetLastError() );
}

TinyMT32::~TinyMT32()
{
	if (d_Random)
		CUDA_SAFE_CALL(cudaFree(d_Random));
}

#include <KMCsplash.h>

void TinyMT32::writeH5(splash::DataCollector* data, int id, const std::string& prefix)
{
	std::string name = prefix + ".tinyMT";

	auto rnd = std::unique_ptr<unsigned int[]>( new unsigned int[m_generatorCount*7] );
	CUDA_SAFE_CALL_THROW( cudaMemcpy(rnd.get(), d_Random, m_generatorCount*sizeof(unsigned int)*7, cudaMemcpyDeviceToHost) );
	data->write(id, KmcSplash::ColTypeUInt32, 2, splash::Selection(splash::Dimensions(m_generatorCount, 7, 1)), name.c_str(), rnd.get());
}

bool TinyMT32::readH5(splash::DataCollector* data, int id, const std::string& prefix)
{
	std::string name = prefix + ".tinyMT";

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
