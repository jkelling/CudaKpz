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

#include "CudaRNG.h"

#include "MersenneTwister_kernel_NVGPUSDK.cu"

#include <cudaError.h>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

#define KMC_L_LCG_N 4
#include <kmcRandom.h>

#define DSFMT_MEXP 19937
#include <dSFMT/dSFMT.h>

extern char _binary_MersenneTwister_dat_start[];

__global__ void RandomGPULCG (unsigned int* m_dRandom, int m_nPerRng);
__global__ void RandomGPULCGS (unsigned int* m_dRandom, int m_nPerRng, unsigned int* seed);

CudaRNG::CudaRNG () : m_dRandom(0), m_nPerRng(0)
{
	memcpy((void*)h_MT, _binary_MersenneTwister_dat_start, sizeof(h_MT));
}

CudaRNG::CudaRNG (const char* file) : m_dRandom(0), m_nPerRng(0)
{
	loadMTGPU(file);
}

CudaRNG::~CudaRNG ()
{
	if(m_dRandom)
		CUDA_SAFE_CALL(cudaFree(m_dRandom));
}

void CudaRNG::setCount (int n)
{
	if(m_dRandom)
		CUDA_SAFE_CALL(cudaFree(m_dRandom));
	m_nPerRng = n/(RNG_COUNT);
	if(n%RNG_COUNT)
		m_nPerRng += 1;
	CUDA_SAFE_CALL(cudaMalloc((void**)&m_dRandom, m_nPerRng*RNG_COUNT*sizeof(unsigned int)));
}

void CudaRNG::generate()
{
	seedMTGPU(rand());
	RandomGPU <<< BLOCKS, THREADS_RED >>> (m_dRandom, m_nPerRng<<L_MT_N_PER_RNG_MUL);
	cudaThreadSynchronize();
	cudaError_t error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		std::cerr << "CUDA error: " << cudaGetErrorString(error) << " in " << __FILE__ << ": " << __LINE__ << '\n';
		exit(1);
	}
}

void CudaRNG::generateLCG()
{
	unsigned int *rnd = new unsigned int[RNG_COUNT];
	if(!rnd)
		throw std::runtime_error("Out of mem in generateLCGS.");
	for(int a = 0; a < RNG_COUNT; ++a)
		rnd[a] = ((double)rand())/RAND_MAX*KMC_LCG_RAND_SUP+1;
	CUDA_SAFE_CALL(cudaMemcpy((void*)m_dRandom, rnd, RNG_COUNT*sizeof(unsigned int), cudaMemcpyHostToDevice));
	delete[] rnd;

	RandomGPULCG <<< BLOCKS, THREADS >>> (m_dRandom, m_nPerRng);
	cudaThreadSynchronize();
	cudaError_t error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		std::cerr << "CUDA error: " << cudaGetErrorString(error) << " in " << __FILE__ << ": " << __LINE__ << '\n';
		exit(1);
	}
}

void CudaRNG::generateLCGS()
{
	unsigned int *rnd = new unsigned int[MT_RNG_COUNT*KMC_LCG_DWORDS];
	if(!rnd)
		throw std::runtime_error("Out of mem in generateLCGS.");
	for(int a = 0; a < MT_RNG_COUNT*KMC_LCG_DWORDS; ++a)
		rnd[a] = (unsigned int)(dsfmt_gv_genrand_close_open()*(KMC_LCG_RAND_SUP))+1;
	unsigned int* seed;
	CUDA_SAFE_CALL(cudaMalloc((void**)&seed, MT_RNG_COUNT*KMC_LCG_DWORDS*sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMemcpy((void*)seed, rnd, MT_RNG_COUNT*KMC_LCG_DWORDS*sizeof(unsigned int), cudaMemcpyHostToDevice));
	delete[] rnd;

	RandomGPULCGS <<< BLOCKS, THREADS_RED >>> (m_dRandom, m_nPerRng<<L_MT_N_PER_RNG_MUL, seed);
	cudaThreadSynchronize();
	cudaError_t error = cudaGetLastError();

	CUDA_SAFE_CALL(cudaFree(seed));
	if(error != cudaSuccess)
	{
		std::cerr << "CUDA error: " << cudaGetErrorString(error) << " in " << __FILE__ << ": " << __LINE__ << '\n';
		exit(1);
	}
}

__global__ void RandomGPULCG (unsigned int* dRandom, int nPerRng)
{
	//seed
	int i = dRandom[blockIdx.x*blockDim.x + threadIdx.x];
	
	const int inc = blockDim.x;
	int gid = blockIdx.x*inc*nPerRng + threadIdx.x;

	for(int a = 0; a < nPerRng; ++a, gid+=inc)
	{
		int y = i/KMC_LCG_Q;
		i = KMC_LCG_A*(i - KMC_LCG_Q*y) - KMC_LCG_R*y;
		if(i<0) i+=KMC_LCG_M;
			dRandom[gid] = i-1;
	}
}

__global__ void RandomGPULCGS (unsigned int* dRandom, int nPerRng, unsigned int* seed)
{
	const int inc = blockDim.x;
	int gid = blockIdx.x*inc*nPerRng + threadIdx.x;
	__shared__ int rng[(CudaRNG::THREADS_RED<<KMC_L_LCG_N)];
	int* j, k, i;
	//seed
	const int BLOCK_RNG_OFFSET = (blockIdx.x * KMC_LCG_DWORDS*blockDim.x); //add block offset
	const int THREAD_RNG_OFFSET = threadIdx.x;
	j = &rng[THREAD_RNG_OFFSET<<KMC_L_LCG_N];
	#pragma unroll
	for (int a = THREAD_RNG_OFFSET; a < (inc<<KMC_L_LCG_N); a+=inc)
		rng[a] = seed[a + BLOCK_RNG_OFFSET];
	k = seed[((KMC_LCG_N)<<inc) + (THREAD_RNG_OFFSET + BLOCK_RNG_OFFSET)];
	i = seed[((KMC_LCG_N+1)<<inc) + (THREAD_RNG_OFFSET + BLOCK_RNG_OFFSET)];

	for(int a = 0; a < nPerRng; ++a, gid+=inc)
	{
		int y = i/KMC_LCG_Q;
		i = KMC_LCG_A*(i - KMC_LCG_Q*y) - KMC_LCG_R*y;
		if(i<0) i+=KMC_LCG_M;
			y = k&KMC_LCG_N_M;
			k = j[y];
			j[y] = i;
			dRandom[gid] = k-1;
	}
}

#include "RC4.cu"

void CudaRNG::generateRC4()
{
	const int RC4SEED_W = BLOCKS*RC4_THREADS*RC4_DWORDS;
	unsigned int *rnd = new unsigned int[RC4SEED_W];
	if(!rnd)
		throw std::runtime_error("Out of mem in generateRC4.");
	for(int a = 0; a < RC4SEED_W; ++a)
	{
		rnd[a] = (unsigned int)(dsfmt_gv_genrand_close_open()*UINT_MAX);
	}
	unsigned int* seed;
	CUDA_SAFE_CALL(cudaMalloc((void**)&seed, RC4SEED_W*sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMemcpy((void*)seed, rnd, RC4SEED_W*sizeof(unsigned int), cudaMemcpyHostToDevice));
	delete[] rnd;

	RandomGPURC4 <<< BLOCKS, RC4_THREADS >>> (m_dRandom, m_nPerRng<<L_RC4_N_PER_RNG_MUL, seed);
	cudaThreadSynchronize();
	cudaError_t error = cudaGetLastError();

	CUDA_SAFE_CALL(cudaFree(seed));
	if(error != cudaSuccess)
	{
		std::cerr << "CUDA error: " << cudaGetErrorString(error) << " in " << __FILE__ << ": " << __LINE__ << '\n';
		exit(1);
	}
}

#include <fstream>
#include "twofish/TwofishKernel.cu"

void CudaRNG::generateTwofish()
{
	const int TWOFISH_SEED_W = BLOCKS*8;
	unsigned int *rnd = new unsigned int[TWOFISH_SEED_W];
	if(!rnd)
		throw std::runtime_error("Out of mem in generateTwofish.");
	for(int a = 0; a < TWOFISH_SEED_W; ++a)
	{
		rnd[a] = (unsigned int)(dsfmt_gv_genrand_close_open()*UINT_MAX);
	}
	unsigned int* seed;
	CUDA_SAFE_CALL(cudaMalloc((void**)&seed, TWOFISH_SEED_W*sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMemcpy((void*)seed, rnd, TWOFISH_SEED_W*sizeof(unsigned int), cudaMemcpyHostToDevice));

	twofishInitConst();

	RandomTwofishECB <<< BLOCKS, TF_THREADS >>> (m_dRandom, m_nPerRng, seed);
	cudaThreadSynchronize();
	cudaError_t error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		std::cerr << "CUDA error: " << cudaGetErrorString(error) << " in " << __FILE__ << ": " << __LINE__ << '\n';
		exit(1);
	}

	delete[] rnd;
	CUDA_SAFE_CALL(cudaFree(seed));
}

bool CudaRNG::twofishTest()
{
	const int BLOCKS = 32;
	const int TWOFISH_SEED_W = BLOCKS*8;
	unsigned int *rnd = new unsigned int[TWOFISH_SEED_W];
	for(int a = 0; a < TWOFISH_SEED_W; ++a)
	{
		rnd[a] = (unsigned int)(dsfmt_gv_genrand_close_open()*UINT_MAX);
	}
	unsigned int* seed;
	CUDA_SAFE_CALL(cudaMalloc((void**)&seed, TWOFISH_SEED_W*sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMemcpy((void*)seed, rnd, TWOFISH_SEED_W*sizeof(unsigned int), cudaMemcpyHostToDevice));

	twofishInitConst();

	const int nPerRng = 1;
	const int data_w = nPerRng*TF_THREADS*BLOCKS;
	unsigned int *data = new unsigned int[data_w];
	for(int a = 0; a < data_w; ++a)
	{
		data[a] = a;//(unsigned int)(dsfmt_gv_genrand_close_open()*UINT_MAX);
	}
	unsigned int* dData;
	CUDA_SAFE_CALL(cudaMalloc((void**)&dData, data_w*sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMemcpy((void*)dData, data, data_w*sizeof(unsigned int), cudaMemcpyHostToDevice));

	encryptTwofishECB <<< BLOCKS, TF_THREADS >>> (dData, nPerRng, seed);
	cudaThreadSynchronize();
	cudaError_t error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		std::cerr << "CUDA error: " << cudaGetErrorString(error) << " in " << __FILE__ << ": " << __LINE__ << '\n';
		exit(1);
	}
	unsigned int *ciphertext = new unsigned int[data_w];
	CUDA_SAFE_CALL(cudaMemcpy(ciphertext, (void*)dData, data_w*sizeof(unsigned int), cudaMemcpyDeviceToHost));
	//for(int a = 0; a < TF_THREADS; ++a)	
	//	std::cerr << data2[a] << '\n';
	decryptTwofishECB <<< BLOCKS, TF_THREADS >>> (dData, nPerRng, seed);
	cudaThreadSynchronize();
	error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		std::cerr << "CUDA error: " << cudaGetErrorString(error) << " in " << __FILE__ << ": " << __LINE__ << '\n';
		exit(1);
	}

	unsigned int *data2 = new unsigned int[data_w];
	CUDA_SAFE_CALL(cudaMemcpy(data2, (void*)dData, data_w*sizeof(unsigned int), cudaMemcpyDeviceToHost));

	int ret = 0;
	for(int a = 0; a < data_w; a+=4)
	{
		int i = a;
		for(; (i < a+4) && (data[i] == data2[i]) ; ++i) ;
		if(i < a+4)
		{
			std::cerr << "error at block " << a/4 << ":\n";
			for(i = a; i< a+4; ++i)
				std::cerr << data[i] << '\t' << data2[i] << "\t(" << ciphertext[i] << ")\n";
			++ret;
			if(ret > 100)
				break;
		}
	}

	CUDA_SAFE_CALL(cudaFree(dData));
	CUDA_SAFE_CALL(cudaFree(seed));
	delete[] data;
	delete[] data2;
	delete[] ciphertext;
	return !ret;
}

void CudaRNG::generateDSFMT()
{
	const int N = RNG_COUNT*m_nPerRng;
	unsigned int *rnd = new unsigned int[N];
	if(!rnd)
		throw std::runtime_error("Out of mem in generateDSFMT.");
	for(int a = 0; a < N; ++a)
		rnd[a] = (unsigned int)(dsfmt_gv_genrand_close_open()*UINT_MAX);
	CUDA_SAFE_CALL(cudaMemcpy((void*)m_dRandom, rnd, N*sizeof(unsigned int), cudaMemcpyHostToDevice));
	delete[] rnd;
}
