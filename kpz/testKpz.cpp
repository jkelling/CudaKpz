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
#include <ctime>
#include <algorithm>

#define DSFMT_MEXP 19937
#include <dSFMT/dSFMT.h>

void mcs(int h[], int N, int n, dsfmt_t& dsfmt, double p, double q)
{
	// std::cerr << p << ' ' << q << " mcs\n";
	for( int a = 0; a < n*N; ++a)
	{
		const int i = dsfmt_genrand_close_open(&dsfmt)*N;
		const double r = dsfmt_genrand_close_open(&dsfmt);
		const int sl = (h[i] - h[(i-1+N)%N]);
		const int sr = (h[(i+1)%N] - h[i]);
		if( sl == -1 && sr == +1)
		{
			if(r < p)
				h[i]+=2;
		}
		else if( sl == +1 && sr == -1)
		{
			if(r < q)
				h[i]-=2;
		}
	}
}

void mcs(int h[], int N, int n, dsfmt_t& dsfmt, const double p[2], const double q[2], bool disorder[])
{
	// std::cerr << p[0] << ' ' << q[0] << " dis\n";
	for( int a = 0; a < n*N; ++a)
	{
		const int i = dsfmt_genrand_close_open(&dsfmt)*N;
		const double r = dsfmt_genrand_close_open(&dsfmt);
		const int sl = (h[i] - h[(i-1+N)%N]);
		const int sr = (h[(i+1)%N] - h[i]);
		// std::cerr << disorder[i] << " diss\n";
		if( sl == -1 && sr == +1)
		{
			if(r < p[disorder[i]])
				h[i]+=2;
		}
		else if( sl == +1 && sr == -1)
		{
			if(r < q[disorder[i]])
				h[i]-=2;
		}
	}
}

std::ostream& autoResponse(int sA[], int sB[], int N, bool disorder[], std::ostream& o, double eps)
{
	double hA = 0, h2A = 0;
	double hB = 0, h2B = 0;
	long long diff = 0;
	for(int a = 0; a < N; ++a)
	{
		hA += sA[a];
		h2A += sA[a]*(double)sA[a];
		hB += sB[a];
		h2B += sB[a]*(double)sB[a];

		if(disorder[a])
			diff -= (sA[a]-sB[a])/eps;
		else
			diff += (sA[a]-sB[a])/eps;
	}
	hA /= N;
	h2A /= N;
	hB /= N;
	h2B /= N;

	const double w2A = h2A - hA*hA;
	const double w2B = h2B - hB*hB;

	const double R = diff/(double)N;
	o << R << '\t' << w2A << '\t' << w2B;
	return o;
}

int main()
{
	static const int N = 1<<13;
	static const double eps = .005;
	static const double p0 = .98;
	dsfmt_t rA, rB;
	dsfmt_init_gen_rand(&rA, time(0));
	int sA[N], sB[N];
	bool disorder[N];

	//init
	for(int a = 0; a < N; a+=2)
	{
		sA[a] = sB[a] = 0;
		sA[a+1] = sB[a+1] = 1;
		disorder[a] = dsfmt_genrand_close_open(&rA)<.5;
		disorder[a+1] = dsfmt_genrand_close_open(&rA)<.5;
	}
	rB = rA;

	static const int waitingTime = 100;
	static const int maxMCS = 1000;
	static const int interval = 10;
	int a = 0;

	const double q0 = 1.-p0;
	const double pd[2] = {p0+eps/2., p0-eps/2.};
	const double qd[2] = {1. - pd[0], 1.-pd[1]};
	while(a < waitingTime)
	{
		const int nextMCS = std::min(interval, waitingTime-a);

		std::cout << a << '\t';
		autoResponse(sA, sB, N, disorder, std::cout, eps) << '\n';
		mcs(sA, N, nextMCS, rA, pd, qd, disorder);
		mcs(sB, N, nextMCS, rB, p0, q0);
		a += nextMCS;
	}
	std::cout << "#waiting time\n";

	std::cout << a << '\t';
	autoResponse(sA, sB, N, disorder, std::cout, eps) << '\n';

	while(a < maxMCS)
	{
		const int nextMCS = std::min(interval, maxMCS-a);

		std::cout << a << '\t';
		autoResponse(sA, sB, N, disorder, std::cout, eps) << '\n';
		mcs(sA, N, nextMCS, rA, p0, q0);
		mcs(sB, N, nextMCS, rB, p0, q0);
		a += nextMCS;
	}

	std::cout << a << '\t';
	autoResponse(sA, sB, N, disorder, std::cout, eps) << '\n';

	for( int a=0; a < N; ++a)
	{
		std::cerr << sA[a] << '\t' << sB[a] << '\n';
	}

	return 0;
}
