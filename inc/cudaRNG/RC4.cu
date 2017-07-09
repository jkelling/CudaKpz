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


const int RC4_DWORDS = 64;
const int RC4_THREADS = 32;

inline static __device__ void swap(char& a, char& b)
{
	const char k = a;
	a = b;
	b = k;
}

__global__ void RandomGPURC4 (unsigned int* dRandom, int nPerRng, unsigned int* seed)
{
	__shared__ char sS[RC4_THREADS<<8];
	char* S = &sS[threadIdx.x<<8];
	for(int a = 0; a < 256; ++a)
		S[a] = a;
	
	const int inc = blockDim.x;
	int gid = blockIdx.x* inc*RC4_DWORDS + threadIdx.x; // offset into seed

	// keysetup
	for(int i = 0, j = 0; i < 256; gid += inc)
	{
		unsigned int t = seed[gid];
		#pragma unroll
		for(int s = 0; s < 4; ++s, ++i, t >>= 8)
		{
			j = (j + S[i] + (t&0xFF))&255;
			swap(S[i], S[j]);
		}
	}

	gid = blockIdx.x* inc*nPerRng + threadIdx.x; // offset into seed
	//keystream (RNG)
	for(int a = 0, i = 0, j = 0; a < nPerRng; ++a, gid += inc)
	{
		unsigned int k = 0;
		#pragma unroll
		for(int s = 0; s < 4; ++s)
		{
			k <<= 8;
			i = (i + 1)&255;
			j = (j + S[i])&255;
			const unsigned char k1 = S[i];
			const unsigned char k2 = S[j];
			S[i] = k2;
			S[j] = k1;
			k |= S[ (k1+k2)&255 ]&255;
		}
		dRandom[gid] = k;
	}
}

