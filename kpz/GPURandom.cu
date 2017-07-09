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

#include "GPURandom.h"

#include "kpzConst.h"
#include <kmcRandom.h>

#include <cassert>

#include <cuda.h>
#include <cudaError.h>

void SLCG64::initCUDA()
{
	CUDA_SAFE_CALL(cudaMalloc((void**)&d_Random, m_generatorCount*sizeof(unsigned long long)));
	if(!d_Random)
	{
		std::cerr << "Insufficient devicememory for random seeds.\n";
		exit(0);
	}
}

void SLCG64::randomize()
{
	assert(m_generatorCount < KPZ_SLCG_SKIP);
	unsigned long long* rnd = new unsigned long long[m_generatorCount];
	rnd[0] = (unsigned long long)(dsfmt_genrand_close_open(m_dsfmt)*((unsigned)-1))
		^ (((unsigned long long)(dsfmt_genrand_close_open(m_dsfmt)*((unsigned)-1)))<<32);
	for (int a = 1; a<m_generatorCount; ++a)
	{
		rnd[a] = SLCGen(rnd[a-1]);
	}
	CUDA_SAFE_CALL( cudaMemcpy(d_Random, rnd, m_generatorCount*sizeof(unsigned long long), cudaMemcpyHostToDevice) );
	delete[] rnd;
}

SLCG64::~SLCG64()
{
	if (d_Random)
		CUDA_SAFE_CALL(cudaFree(d_Random));
}
