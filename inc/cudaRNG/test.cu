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

#include <iostream>
#include <cuda.h>
#include <cstdlib>
#include <ctime>
#include <cstring>

#include "CudaRNG.h"
#include <fft.h>

#define DSFMT_MEXP 19937
#include <dSFMT/dSFMT.h>

void corr (unsigned int* d, int N, std::ostream& o)
{
	complex c[N];
	for(int a = 0; a < N; ++a)
	{
		c[a].Re = (double)d[a]/UINT_MAX;
		c[a].Im = 0;
	}
	/*
	fft(c, N, 1);
	norm(c, N, sqrt(N));
	absSq(c, N);
	fft(c, N, -1);
	norm(c, N, sqrt(N));*/
	for(int a = 0; a < N; ++a)
		o << a << '\t' << c[a].Re << '\t' << c[a].Im << '\n';
	/*
	if(maxDelta == -1 || maxDelta >= N)
		maxDelta = N-1;
	for(int delta = 1; delta < maxDelta; ++delta)
	{
		double sum;
		for(int a = 0; a < N; ++a)
		{
			const register double diff = d[a] - d[(a+delta)%N];
			sum += diff*diff;
		}
		sum /= N;
		o << delta << '\t' << 1/sum << '\n';
	}*/
}

std::ostream& put(unsigned int* d, int N, int l = 10)
{
	for(int a = 0; a < N; ++a)
	{
		std::cout << d[a] << ' ';
		if(!((a+1)%l))
			std::cout << '\n';
	}
	return std::cout;
}

int main(int argc, char* argv[])
{
	srand(time(0));
	dsfmt_gv_init_gen_rand(rand());
	for( int a = 0; a < 100000; ++a)
	{
		dsfmt_gv_genrand_close_open();
	}

	int N = 4096;
	if(argc > 1)
		N = atoi(argv[1]);
	std::cerr << N << '\n';

	CudaRNG rnd;
	rnd.setCount(N);
	N = rnd.n();
	std::cerr << N << '\n';
	unsigned int a[N];

	//rnd.twofishTest();
//	for(int a = 0; a < 400; ++a)
		rnd.generateLCG();
	cudaMemcpy((void*)a, rnd.dRandom(), N*sizeof(unsigned int), cudaMemcpyDeviceToHost); 
	corr(a, N, std::cout);
	//rnd.generateLCGS();
	//cudaMemcpy((void*)a, rnd.dRandom(), N*sizeof(unsigned int), cudaMemcpyDeviceToHost); 

	return 0;
}
