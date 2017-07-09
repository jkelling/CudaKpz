typedef unsigned int u32;
typedef unsigned char BYTE;

__constant__ BYTE RS[4][8];
__constant__ BYTE Q0[256];
__constant__ BYTE Q1[256];

#include "tables.h"
#include <cudaError.h>
#include <cuda.h>

static void twofishInitConst()
{
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(RS, (void*)RS_HOST, sizeof(RS_HOST)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(Q0, (void*)Q0_HOST, sizeof(Q0_HOST)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(Q1, (void*)Q1_HOST, sizeof(Q1_HOST)));
}

//#define TWOFISH_SMALL_BLOCK // allows blockDim < 256, still >= 20
#define TWOFISH_MULTI // turn off let on thread do keysetup alone
#include "TwofishInternal.cu"

const int TF_BLOCK_SIZE_W = 4;
const int TF_THREADS = 512;

/* Always call kernel with a minimum blockDim of 256, fullKey() breaks otherwise. */
__global__ void encryptTwofishECB (u32* dData, int nPerRng, u32* dKey)
{
	__shared__ u32 K[40];
	__shared__ u32 QF[4][256];
	__shared__ u32 data[TF_THREADS*TF_BLOCK_SIZE_W];
	
	init(dKey, K, QF);
	__syncthreads();
	
	const int inc = blockDim.x*TF_BLOCK_SIZE_W;
	int gid = blockIdx.x * nPerRng*inc;
	const int lid = threadIdx.x*TF_BLOCK_SIZE_W;

	for(int a = 0; a < nPerRng; ++a, gid += inc)
	{
		#pragma unroll
		for(int i = threadIdx.x; i < inc; i += blockDim.x)
			data[i] = dData[gid + i];
		__syncthreads();
		
		encryptblock(K, QF, &(data[lid]));
		
		__syncthreads();
		#pragma unroll
		for(int i = threadIdx.x; i < inc; i += blockDim.x)
			dData[gid + i] = data[i];
		__syncthreads();
	}
}

__global__ void decryptTwofishECB (u32* dData, int nPerRng, u32* dKey)
{
	__shared__ u32 K[40];
	__shared__ u32 QF[4][256];
	__shared__ u32 data[TF_THREADS*TF_BLOCK_SIZE_W];
	
	init(dKey, K, QF);
	__syncthreads();
	
	const int inc = blockDim.x*TF_BLOCK_SIZE_W;
	int gid = blockIdx.x * nPerRng*inc;
	const int lid = threadIdx.x*TF_BLOCK_SIZE_W;
	
	for(int a = 0; a < nPerRng; ++a, gid += inc)
	{
		#pragma unroll
		for(int i = threadIdx.x; i < inc; i += blockDim.x)
			data[i] = dData[gid + i];
		__syncthreads();
		
		decryptblock(K, QF, &(data[lid]));
		
		__syncthreads();
		#pragma unroll
		for(int i = threadIdx.x; i < inc; i += blockDim.x)
			dData[gid + i] = data[i];
		__syncthreads();
	}
}

__global__ void RandomTwofishECB (u32* dRandom, int nPerRng, u32* dKey)
{
	__shared__ u32 K[40];
	__shared__ u32 QF[4][256];
	__shared__ u32 data[TF_THREADS*TF_BLOCK_SIZE_W];
	
	init(dKey, K, QF);
	__syncthreads();
	
	const int inc = blockDim.x*TF_BLOCK_SIZE_W;
	int gid = blockIdx.x * (nPerRng>>2)*inc;
	const int lid = threadIdx.x*TF_BLOCK_SIZE_W;
	int counter = 0;

	for(int a = 0; a < (nPerRng>>2); ++a, gid += inc, ++counter)
	{
		data[lid] = counter;
		data[lid + 1] = threadIdx.x;
		data[lid + 2] = blockIdx.x;
		data[lid + 3] = 0;
		__syncthreads();
		
		encryptblock(K, QF, &(data[lid]));
		
		__syncthreads();
		#pragma unroll
		for(int i = threadIdx.x; i < inc; i += blockDim.x)
			dRandom[gid + i] = data[i];
		__syncthreads();
	}
}

