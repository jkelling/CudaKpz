/***************************************************************************
*   Copyright 2011 - 2012 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef KPZ_KPZ_COPY
#define KPZ_KPZ_COPY

//#define SINGLE_THREAD_COPY
#define KMC_BULK_CPY

#ifndef __OPENCL_VERSION__
template<class LocalLayout>
static inline __device__ void copyIn (unsigned int* system, unsigned int* d_System, const int blocks, int xw, int yw/*, int* comm*/)
#else
inline void copyIn(__local unsigned int* system, __global unsigned int* d_System, const int blocks, int xw, int yw)
#endif
{
	#ifdef KMC_BULK_CPY
		#define cpy(a,b) \
			system[ a ] = d_System[ b ];
		int sindex = GPU_THREAD_ID_X + (GPU_THREAD_ID_Y<<(LocalLayout::L_BLOCK_DIM_X_W+LocalLayout::L_THREAD_CELL_DIM_Y_W));
		xw = ((xw + (blocks << LocalLayout::L_BLOCK_DIM_X_W) + GPU_THREAD_ID_X) & SIZE_M_DIM_X_W);
		yw = ((yw + ((blocks>>SIZE_L_BLOCKS_X) << LocalLayout::L_BLOCK_DIM_Y_W) + (GPU_THREAD_ID_Y<<LocalLayout::L_THREAD_CELL_DIM_Y_W)) & SIZE_M_DIM_Y_W);
		for(int a = 0; a < LocalLayout::THREAD_CELL_DIM_Y_W; ++a)
		{
			int gx = xw;
			#pragma unroll
			for(int x = 0; x < LocalLayout::THREAD_CELL_DIM_X_W; ++x)
			{
				//cpy(sindex, gx + (yw<<SIZE_L_DIM_X_W));
				system[sindex] = d_System[gx + (yw<<SIZE_L_DIM_X_W)];
				gx = ((gx + LocalLayout::THREADS_X) & SIZE_M_DIM_X_W);
				sindex += LocalLayout::THREADS_X;
			}
			yw = ((yw + 1) & SIZE_M_DIM_Y_W);
		}
		#undef cpy
	#else
		xw = ((xw + (blocks << LocalLayout::L_BLOCK_DIM_X_W)) & SIZE_M_DIM_X_W);
		yw = ((yw + ((blocks>>SIZE_L_BLOCKS_X) << LocalLayout::L_BLOCK_DIM_Y_W)) & SIZE_M_DIM_Y_W);
		#ifdef SINGLE_THREAD_COPY
		if(GPU_THREAD_ID_X | GPU_THREAD_ID_Y)
			return;
		#else
		if(GPU_THREAD_ID_Y > 1) // for some reason this fixes the inconsistency arising for large numbers of threads.
			return;
		#endif
		for (int a = 0, y = yw; a<LocalLayout::BLOCK_DIM_Y_W; ++a, (++y)&=SIZE_M_DIM_Y_W)
		{
			#ifdef SINGLE_THREAD_COPY
			for (int s = 0, x = (xw)&SIZE_M_DIM_X_W
				; s < LocalLayout::BLOCK_DIM_X_W; s+=1, (x+=1)&=SIZE_M_DIM_X_W)
			#else
			for (int s = (GPU_THREAD_ID_Y<<LocalLayout::L_THREADS_X)+GPU_THREAD_ID_X
				, x = (xw+((GPU_THREAD_ID_Y<<LocalLayout::L_THREADS_X)+GPU_THREAD_ID_X))&SIZE_M_DIM_X_W
				; s < LocalLayout::BLOCK_DIM_X_W
				; s+=LocalLayout::THREADS, (x+=LocalLayout::THREADS)&=SIZE_M_DIM_X_W)
			#endif
			{
				system[(a<<LocalLayout::L_BLOCK_DIM_X_W) + s] = d_System[(y<<SIZE_L_DIM_X_W) + x];
			}
			#ifndef SINGLE_THREAD_COPY
			//__syncthreads();
			#endif
		}
	#endif
}

#ifndef __OPENCL_VERSION__
template<class LocalLayout>
static inline __device__ void copyOut (unsigned int* system, unsigned int* d_System, const int blocks, int xw, int yw)
#else
inline void copyOut(__local unsigned int* system, __global unsigned int* d_System, const int blocks, int xw, int yw)
#endif
{
	#ifdef KMC_BULK_CPY
		#define cpy(a,b) \
			d_System[ b ] = system[ a ];
		int sindex = GPU_THREAD_ID_X + (GPU_THREAD_ID_Y<<(LocalLayout::L_BLOCK_DIM_X_W+LocalLayout::L_THREAD_CELL_DIM_Y_W));
		xw = ((xw + (blocks << LocalLayout::L_BLOCK_DIM_X_W) + GPU_THREAD_ID_X) & SIZE_M_DIM_X_W);
		yw = ((yw + ((blocks>>SIZE_L_BLOCKS_X) << LocalLayout::L_BLOCK_DIM_Y_W) + (GPU_THREAD_ID_Y<<LocalLayout::L_THREAD_CELL_DIM_Y_W)) & SIZE_M_DIM_Y_W);
		for(int a = 0; a < LocalLayout::THREAD_CELL_DIM_Y_W; ++a)
		{
			int gx = xw;
			#pragma unroll
			for(int x = 0; x < LocalLayout::THREAD_CELL_DIM_X_W; ++x)
			{
				//cpy(sindex, gx + (yw<<KPZ_DEVICE_SYSTEM_SIZE.m_lDimXW));
				d_System[gx + (yw<<SIZE_L_DIM_X_W)] = system[sindex];
				gx = ((gx + LocalLayout::THREADS_X) & SIZE_M_DIM_X_W);
				sindex += LocalLayout::THREADS_X;
			}
			yw = ((yw + 1) & SIZE_M_DIM_Y_W);
		}
		#undef cpy
	#else
		xw = ((xw + (blocks << LocalLayout::L_BLOCK_DIM_X_W)) & SIZE_M_DIM_X_W);
		yw = ((yw + ((blocks>>SIZE_L_BLOCKS_X) << LocalLayout::L_BLOCK_DIM_Y_W)) & SIZE_M_DIM_Y_W);
		#ifdef SINGLE_THREAD_COPY
		if(GPU_THREAD_ID_X | GPU_THREAD_ID_Y)
			return;
		#else
		if(GPU_THREAD_ID_Y > 1) // for some reason this fixes the inconsistency arising for large numbers of threads.
			return;
		#endif
		for (int a = 0, y = yw; a<LocalLayout::BLOCK_DIM_Y_W; ++a, (++y)&=SIZE_M_DIM_Y_W)
		{
			#ifdef SINGLE_THREAD_COPY
			for (int s = 0, x = (xw)&SIZE_M_DIM_X_W
				; s < LocalLayout::BLOCK_DIM_X_W; s+=1, (x+=1)&=SIZE_M_DIM_X_W)
			#else
			for (int s = (GPU_THREAD_ID_Y<<LocalLayout::L_THREADS_X)+GPU_THREAD_ID_X
				, x = (xw+((GPU_THREAD_ID_Y<<LocalLayout::L_THREADS_X)+GPU_THREAD_ID_X))&SIZE_M_DIM_X_W
				; s < LocalLayout::BLOCK_DIM_X_W
				; s+=LocalLayout::THREADS, (x+=LocalLayout::THREADS)&=SIZE_M_DIM_X_W)
			#endif
			{
				d_System[(y<<SIZE_L_DIM_X_W) + x] = system[(a<<LocalLayout::L_BLOCK_DIM_X_W) + s];
			}
			#ifndef SINGLE_THREAD_COPY
			//__syncthreads();
			#endif
		}
	#endif
}

template<class Cpy, class LocalLayout>
void __device__ simpleCopyDT(const Cpy& cpy, const int& blocks, const int& block, const int xw, const int yw)
{
	const int Y = (blocks>>KPZ_DEVICE_SYSTEM_SIZE.m_lBlocksX<<(LocalLayout::L_BLOCK_DIM_Y_W+1))+yw
		+ ((block&2) ? LocalLayout::BLOCK_DIM_Y_W : 0);
	const int X = ((blocks&((1<<KPZ_DEVICE_SYSTEM_SIZE.m_lBlocksX)-1))<<(LocalLayout::L_BLOCK_DIM_X_W+1))+xw
		+ ((block&1) ? LocalLayout::BLOCK_DIM_X_W : 0);
	for(int i = threadIdx.x + (threadIdx.y<<LocalLayout::L_THREADS_X); i < LocalLayout::BLOCK_DATA; i += LocalLayout::THREADS)
	{
		int li = i;
		const int x = (li % (LocalLayout::BLOCK_DIM_X_W+1) + X)&KPZ_DEVICE_SYSTEM_SIZE.m_mDimXW;
		li /= (LocalLayout::BLOCK_DIM_X_W+1);
		const int y = (li + Y)&KPZ_DEVICE_SYSTEM_SIZE.m_mDimYW;
		const int gi = x + (y<<KPZ_DEVICE_SYSTEM_SIZE.m_lDimXW);
		cpy(i, gi);
	}
}

struct CopyIn {
	unsigned int* d_System, *system;
	inline void __device__ operator() (const int a, const int b) const {
		system[ a ] = d_System[ b ];
	}
};
struct CopyOut {
	unsigned int* d_System, *system;
	inline void __device__ operator() (const int a, const int b) const {
		d_System[ b ] = system[ a ];
	}
};

template<class LocalLayout>
static inline __device__ void copyInDT (unsigned int* system, unsigned int* d_System, const int blocks, int xw, int yw, int block)
{
	CopyIn cpy = {d_System, system};
	simpleCopyDT<CopyIn, LocalLayout>(cpy, blocks, block, xw, yw);
}

template<class LocalLayout>
static inline __device__ void copyOutDT (unsigned int* system, unsigned int* d_System, const int blocks, int xw, int yw, int block)
{
	CopyOut cpy = {d_System, system};
	simpleCopyDT<CopyOut, LocalLayout>(cpy, blocks, block, xw, yw);
}

#ifdef SINGLE_THREAD_COPY
	#undef SINGLE_THREAD_COPY
#endif
#endif
