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

#include "Prand.h"

#include <cuda.h>
#include <cudaError.h>

#include <iostream>

void PrandMT::initCUDA()
{
	mt19937_init_device_consts_();

	size_t free, total;
	cuMemGetInfo(&free, &total);
	if(free < 64 || free < m_generatorCount)
	{
		std::cerr << "Insufficient devicememory for random numbers.\n";
		exit(0);
	}
	m_bufferSize = 1<<PRNG_L_MAX_BUFFER_SIZE;
	while(free < m_bufferSize)
		m_bufferSize >>= 1;
	CUDA_SAFE_CALL(cudaMalloc((void**)&d_Random, m_bufferSize));
	m_bufferSize /= sizeof(unsigned int);
}

PrandMT::PrandMT(size_t generatorCount, dsfmt_t* dsfmt)
	: m_generatorCount(generatorCount), m_dsfmt(dsfmt)
{
	m_state = new mt19937_state;
	initCUDA();
	m_used = m_bufferSize;
}

PrandMT::~PrandMT()
{
	delete m_state;
	mt19937_free_device_consts_();
	if (d_Random)
		CUDA_SAFE_CALL(cudaFree(d_Random));
}

void PrandMT::randomize()
{
	mt19937_init_sequence_(m_state, (unsigned long long)(dsfmt_genrand_close_open(m_dsfmt)*((unsigned)-1))); 
	fillBuffer(0);
}

bool PrandMT::fillBuffer(int numbers)
{
	if(numbers > m_bufferSize)
	{
		std::cerr << "Cannot generate as many random number as requested in one run. Increase L_MCS_DIV.\n";
		return false;
	}
	mt19937_generate_gpu_array_(m_state,d_Random,m_bufferSize);
	mt19937_skipahead_(m_state, m_bufferSize, 1);
	m_used = 0;
	return true;
}
