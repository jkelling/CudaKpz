/***************************************************************************
*   Copyright 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#include "correlator.h"
#include "kpzConst.h"

#include <kmcExceptCUDA.h>

#include <cuda.h>
#include <cudaError.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>

const int MAX_THREADS = 1024;
const int N_META = 2;
const int N_RESULTS = Kpz::Roughness::H_ORDERS + 2;

/* dimXW must be multiple of blockDim.x
 * pmeta[0] = hSys
 * pmeta[1] = hSnap
 */
__global__ void correlateKernelCBBit(const unsigned int* d_sys, const unsigned int* d_snap, const int* pmeta
		, const unsigned int dimXW, const unsigned int quadrant, const unsigned int ymin, int* d_results
		// , int* d_heightmap
		);

Kpz::Roughness Kpz::Correlator::correlateCuda(unsigned int* d_System, int nBlocks)
{
	if(m_scheduler->encoding() != SchedulerService::ENC_CBBIT)
		return correlateThreaded();

	const SystemSize& size = m_scheduler->size();
	const unsigned int* system = m_scheduler->system();
	const int threads = std::min(MAX_THREADS, size.dimX()>>6);
	// std::cerr << "correlateCuda: running " << nBlocks << " blocks and " << threads << " threads per block\n";

	static const int N_CHUNKS = 4;
	// init all device stuff
	cudaError_t cudaError = cudaHostRegister(m_snapshot, size.sizeW()<<2, cudaHostRegisterPortable);
	bool snapshotPinned = true;
	if(cudaError != cudaSuccess)
	{
		snapshotPinned = false;
		const char* string = cudaGetErrorString(cudaError);
		std::cerr << "[Correlator::correlateCuda][WW] Error when trying to pin host memory for snapshot:\n\t" << string << '\n';
		cudaGetLastError();
	}
	int* d_results;
	CUDA_SAFE_CALL(cudaMalloc((void**)&d_results, size.dimY()*(N_RESULTS*sizeof(double))));
	unsigned int* d_snap;
	const size_t snapMemPerChunkW = (size.dimX()>>4)*nBlocks;
	CUDA_SAFE_CALL(cudaMalloc((void**)&d_snap, (snapMemPerChunkW*N_CHUNKS*sizeof(unsigned int))));
	const size_t metaPerChunkW = N_META*nBlocks;
	int* d_meta, *meta = new int[N_CHUNKS*metaPerChunkW];
	CUDA_SAFE_CALL(cudaMalloc((void**)&d_meta, (metaPerChunkW*N_CHUNKS*sizeof(int))));
	CUDA_SAFE_CALL(cudaHostRegister(meta, (metaPerChunkW*N_CHUNKS*sizeof(int)), cudaHostRegisterPortable));
	cudaStream_t streams[N_CHUNKS];
	cudaEvent_t memTransferEvents[N_CHUNKS];
	for(int a = 0; a < N_CHUNKS; ++a)
	{
		CUDA_SAFE_CALL(cudaStreamCreate(&streams[a]));
		CUDA_SAFE_CALL(cudaEventCreate(&memTransferEvents[a]));
	}

	int hSys = 0, hSnap = 0; // init with hight left of (0,-1)
	{
		auto s = size.slopeX_CBBit(0, size.dimY()-1);
		hSys = (s(system)) ? -1 : +1; // go back/left one site in x-dir
		hSnap = (s(m_snapshot)) ? -1 : +1;
	}

		// int* d_heightmap;
		// const size_t memHeightmapW = size.size();
		// CUDA_SAFE_CALL(cudaMalloc((void**)&d_heightmap, (memHeightmapW*sizeof(int))));
		// int* heightmap = new int[memHeightmapW];
		// for(int a = 0; a < memHeightmapW; ++a)
		// 	heightmap[a] = INT_MAX;
		// CUDA_SAFE_CALL(cudaMemcpy(d_heightmap, heightmap, (memHeightmapW)*sizeof(int), cudaMemcpyHostToDevice));

	int stream = 0;
	for(unsigned int ymin = 0; ymin < size.dimY(); ymin +=nBlocks)
	{
		int blocks = std::min(ymin+nBlocks, (unsigned)size.dimY()) - ymin;

		{
			const size_t memPerChunkQuadrant = (size.dimX()>>4)*blocks;
			for(int a = 0; a < 4; ++a)
			{
				CUDA_SAFE_CALL(cudaMemcpyAsync(
					(void*)(d_snap + stream*snapMemPerChunkW + a*(memPerChunkQuadrant>>2))
					, (void*)&m_snapshot[(ymin<<(size.lDimX()-6)) + a*(size.sizeW()>>2)], memPerChunkQuadrant
					, cudaMemcpyHostToDevice, streams[stream]));
			}
		}

		int* localMeta = &meta[stream*metaPerChunkW];
		const size_t indexBase = (size.sizeW()>>2) + (size.dimX()>>6)-1;
		cudaEventSynchronize(memTransferEvents[stream]);
		for(int y = ymin; y < ymin+blocks; ++y)
		{
			const size_t index = indexBase + (y<<(size.lDimX()-6)) + ((y&1) ? 0 : (size.sizeW()>>1));
			hSys += ((system[index]>>31)&1) ? +1 : -1;
			hSnap += ((m_snapshot[index]>>31)&1) ? +1 : -1;
			localMeta[((y-ymin)<<1)] = hSys;
			localMeta[((y-ymin)<<1) + 1] = hSnap;
			// std::cerr << "Fe\t" << size.dimX()-1 << '\t' << y << "\t" << localMeta[((y-ymin)<<1)+1] << '\t' << hSys<<'\n';
		}
		CUDA_SAFE_CALL(cudaMemcpyAsync(d_meta + metaPerChunkW*stream, localMeta, metaPerChunkW*sizeof(int)
			, cudaMemcpyHostToDevice, streams[stream]));
		cudaEventRecord(memTransferEvents[stream]);
		correlateKernelCBBit<<< blocks, threads, 0, streams[stream]>>>
			(d_System, d_snap + (stream*snapMemPerChunkW), d_meta+metaPerChunkW*stream, size.dimX()>>6, size.sizeW()>>2
			 , ymin, d_results + (ymin*(N_RESULTS*2))
			 // , d_heightmap
			 );

		stream = (stream+1)%N_CHUNKS;
	}

	for(int a = 0; a < N_CHUNKS; ++a)
	{
		CUDA_SAFE_CALL(cudaStreamSynchronize(streams[a]));
		CUDA_SAFE_CALL(cudaStreamDestroy(streams[a]));
		CUDA_SAFE_CALL(cudaEventDestroy(memTransferEvents[a]));
	}

	// {
	// 	CUDA_SAFE_CALL(cudaMemcpy(heightmap, d_heightmap, (memHeightmapW)*sizeof(int), cudaMemcpyDeviceToHost));
	// 	CUDA_SAFE_CALL(cudaFree(d_heightmap));
	// 	std::ostringstream os;
	// 	static int mcs = 0;
	// 	++mcs;
	// 	os << "gpuhm_" << mcs << ".xyz";
	// 	std::ofstream of(os.str().c_str());
	// 	for(int y = 0; y < size.dimY(); ++y)
	// 	{
	// 		for(int x = 1; x < size.dimX(); ++x)
	// 		{
	// 			const auto h = heightmap[y*size.dimX()+x];
	// 			of << "Fe\t" << x << '\t' << y << '\t' << h << '\t' << h << '\n';
	// 		}
	// 		const auto h = heightmap[y*size.dimX()];
	// 		of << "Fe\t" << 0 << '\t' << y << '\t' << h << '\t' << h << '\n';
	// 	}
	// }


	double* results = new double[size.dimY()*N_RESULTS];
	CUDA_SAFE_CALL(cudaMemcpy(results, d_results, size.dimY()*(N_RESULTS*sizeof(double)), cudaMemcpyDeviceToHost));

	Roughness r;
	m_ch = m_cs = 0.;
	for(int y = 0; y < size.dimY(); ++y)
	{
		double *res = results + y*N_RESULTS;
		for(int a = 0; a < Kpz::Roughness::H_ORDERS; ++a)
			r.h[a] += res[a];
		m_ch += res[Kpz::Roughness::H_ORDERS];
		m_cs += res[Kpz::Roughness::H_ORDERS+1];
	}
	r.normAndSetW2(size.size());
	m_ch = m_ch/size.size() - r.h[0]*m_h;
	m_cs = m_cs/size.size()/2.;

	delete[] results;

	if(snapshotPinned)
		CUDA_SAFE_CALL(cudaHostUnregister(m_snapshot));
	CUDA_SAFE_CALL(cudaHostUnregister(meta));
	CUDA_SAFE_CALL(cudaFree(d_snap));
	CUDA_SAFE_CALL(cudaFree(d_results));
	CUDA_SAFE_CALL(cudaFree(d_meta));
	delete[] meta;

	cudaError = cudaGetLastError();
	if(cudaError != cudaSuccess)
	{
		throw KmcExcept::CUDAError(cudaError, __FILE__, __LINE__);
	}

	return r;
}

__device__ inline void ffHeightCBBit(int& h, unsigned int* tx)
{
	for(int x = 0; x < 32; ++x)
	{
		 h += ((tx[0]>>x)&1) ? +1 : -1;
		 h += ((tx[1]>>x)&1) ? +1 : -1;
	}
}

__device__ inline void sumXnor(int& tcs, const unsigned int data)
{
	for(int x = 0; x < 32; ++x)
		tcs += ((data>>x)&1) ? -1 : +1;
}

union SplitDouble {
	double val;
	struct {
		int lo, hi;
	};

	__device__ inline SplitDouble(const double& val = 0.) : val(val) {}
	__device__ inline SplitDouble(const int lo, const int hi) : lo(lo), hi(hi) {}
};

/* blockDim.x must be power of two */
__device__ inline void reduceDoubleFermi(SplitDouble& val, int* sharedBuffer)
{
	__syncthreads();
	sharedBuffer[threadIdx.x] = val.lo;
	sharedBuffer[threadIdx.x + MAX_THREADS] = val.hi;
	for(int n = blockDim.x>>1; n > 0; n>>=1)
	{
		__syncthreads();
		if(threadIdx.x < n)
		{
			SplitDouble tmp;
			tmp.lo = sharedBuffer[threadIdx.x+n];
			tmp.hi = sharedBuffer[(threadIdx.x+n) + MAX_THREADS];
			val.val += tmp.val;
			sharedBuffer[threadIdx.x] = val.lo;
			sharedBuffer[threadIdx.x + MAX_THREADS] = val.hi;
		}
	}
}

/* blockDim.x must be power of two */
__device__ inline void reduceIntFermi(int& val, int* sharedBuffer)
{
	__syncthreads();
	sharedBuffer[threadIdx.x] = val;
	for(int n = blockDim.x>>1; n > 0; n>>=1)
	{
		__syncthreads();
		if(threadIdx.x < n)
		{
			val += sharedBuffer[threadIdx.x+n];
			sharedBuffer[threadIdx.x] = val;
		}
	}
}

__device__ inline void addH(double hSys, SplitDouble* h)
{
	h[0].val += hSys;
	double acc = hSys;
#pragma unroll
	for(int a = 1; a < Kpz::Roughness::H_ORDERS; ++a)
		h[a].val += (acc*=hSys);
}

__global__ void correlateKernelCBBit(const unsigned int* d_sys, const unsigned int* d_snap, const int* pmeta
		, const unsigned int dimXW, const unsigned int quadrant, unsigned int ymin, int* d_results
		// , int* d_heightmap
		)
{
	ymin += blockIdx.x;
	__shared__ int meta[N_META];
	__shared__ int rightHeight[MAX_THREADS*2+1]; // pad to eliminate bank conflicts
	if(threadIdx.x < N_META)
		meta[threadIdx.x] = pmeta[blockIdx.x*N_META + threadIdx.x];

	/*SplitDouble h = 0., h2 = 0., ch = 0.;*/
	SplitDouble h[Kpz::Roughness::H_ORDERS], ch = 0.;
#pragma unroll
	for(int a = 0; a < Kpz::Roughness::H_ORDERS; ++a)
		h[a].val = 0.;
	int tcs = 0;
	for(int x = threadIdx.x; x < dimXW; x+= blockDim.x)
	{
		unsigned int sysTX[2], snapTX[2];
		sysTX[0] = d_sys[x + dimXW*(ymin)];
		sysTX[1] = d_sys[x + dimXW*(ymin) + (quadrant<<1)];
		int hSys = 0;
		ffHeightCBBit(hSys, sysTX);
		rightHeight[threadIdx.x] = hSys;
		snapTX[0] = d_snap[x + dimXW*blockIdx.x]; // even x
		snapTX[1] = d_snap[x + dimXW*(blockIdx.x + (gridDim.x<<1))]; // odd x
		int hSnap = 0;
		ffHeightCBBit(hSnap, snapTX);
		rightHeight[MAX_THREADS + 1 + threadIdx.x] = hSnap;
		__syncthreads();

		if(threadIdx.x < 2)
		{
			hSys = meta[threadIdx.x]; // get initial height of this line
			int* ptr = rightHeight + ((threadIdx.x) ? MAX_THREADS +1 : 0);
			for(int a = 0; a < blockDim.x; ++a)
			{
				hSnap = ptr[a];
				ptr[a] = hSys;
				hSys += hSnap;
			}
			meta[threadIdx.x] = hSys; // for the next run
		}
		__syncthreads();

		hSys = rightHeight[threadIdx.x];
		hSnap = rightHeight[MAX_THREADS+1+threadIdx.x];
		for(int xb = 0; xb < 32; ++xb)
		{
			hSys += ((sysTX[ymin&1]>>xb)&1) ? +1 : -1;
			hSnap += ((snapTX[ymin&1]>>xb)&1) ? +1 : -1;
			// d_heightmap[ymin*dimXW*64 + x*64+xb*2] = hSys;
			addH(hSys, h);
			ch.val += hSnap*(double)hSys;
			hSys += ((sysTX[(~ymin)&1]>>xb)&1) ? +1 : -1;
			hSnap += ((snapTX[(~ymin)&1]>>xb)&1) ? +1 : -1;
			// d_heightmap[ymin*dimXW*64 + x*64+xb*2+1] = hSys;
			addH(hSys, h);
			ch.val += hSnap*(double)hSys;
		}

#if 1
		sysTX[0] ^= snapTX[0];
		sysTX[1] ^= snapTX[1];
		snapTX[0] = d_snap[x + dimXW*(blockIdx.x + gridDim.x)]; // even y
		snapTX[1] = d_snap[x + dimXW*(blockIdx.x + (gridDim.x<<1) + gridDim.x)]; // odd y
		sumXnor(tcs, sysTX[0]);
		sumXnor(tcs, sysTX[1]);

		sysTX[0] = d_sys[x + dimXW*ymin + quadrant]; // even y
		sysTX[1] = d_sys[x + dimXW*ymin + (quadrant<<1) + quadrant]; // odd y
		sysTX[0] ^= snapTX[0];
		sysTX[1] ^= snapTX[1];
		sumXnor(tcs, sysTX[0]);
		sumXnor(tcs, sysTX[1]);
#endif
	}

	for(int a = 0; a < Kpz::Roughness::H_ORDERS; ++a)
		reduceDoubleFermi(h[a], rightHeight);
	reduceDoubleFermi(ch, rightHeight);
	reduceIntFermi(tcs, rightHeight);
	if(threadIdx.x == 0)
	{
		for(int a = 0; a < Kpz::Roughness::H_ORDERS; ++a)
			((double*)rightHeight)[a] = h[a].val;
		((double*)rightHeight)[Kpz::Roughness::H_ORDERS] = ch.val;
		((double*)rightHeight)[Kpz::Roughness::H_ORDERS+1] = (double)tcs;
	}
	__syncthreads();
	if(threadIdx.x < N_RESULTS*2)
		d_results[(blockIdx.x*N_RESULTS*2) + threadIdx.x] = rightHeight[threadIdx.x];
}

#if 0
__global__ void correlateKernelCO(const unsigned int* d_sys, const unsigned int* d_snap, const int* pmeta, const int dimXW)
{
	__shared__ int meta[N_META];
	__shared__ int rightHeightSys[MAX_THREADS*2+1]; // pad to avoid bank conflicts
	int* rightHeightSnap = rightHeightSys+MAX_THREADS+1;
	if(threadIdx.x < N_META)
		meta[threadIdx.x] = pmeta[blockIdx.x*N_META + threadIdx.x];
	unsigned int* sys, snap;
	for(int x = threadIdx.x; x < dimXW; x+= blockDim.x)
	{
		sys = s_sys[dimXW*blockIdx.x + x];
		int hSys = 0;
		for(int x = 0; x < BASIC_CELL_DIM_X; ++x)
			 hSys += ((sys>>(x<<1))&1) ? +1 : -1;
		rightHeight[threadIdx.x] = hSys;
		snap = s_snap[dimXW*blockIdx.x + x];
		int hSnap = 0;
		for(int x = 0; x < BASIC_CELL_DIM_X; ++x)
			 hSnap += ((snap>>(x<<1))&1) ? +1 : -1;
		rightHeight[MAX_THREADS + 1 + threadIdx.x] = hSnap;
		__syncthreads();

		if(threadIdx.x < 2)
		{
			int runningSum = meta[threadIdx.x]; // get initial height of this line
			int* ptr = (threadIdx.x) ? rightHeightSnap : rightHeightSys;
			for(int a = 0; a < blockDim.x; ++a)
				runningSum = (ptr[a]+=runningSum);
			meta[threadIdx.x] = runningSum; // for the next run
		}
		__syncthreads();

		hSys = rightHeightSys[threadIdx.x];
		hSnap = rightHeightSnap[threadIdx.x];
		int shift = 0;
		for(int y = 0; 
	}
}
#endif
