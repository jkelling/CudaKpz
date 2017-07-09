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

#define DSFMT_MEXP 19937
#include <dSFMT/dSFMT.h>

#include <benchmarking.h>

#include "lattice.h"
#include "dim.h"
#include "sumSim.h"
#include "slopeSim.h"

#include <iostream>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <ctime>

#ifndef DIMENSION
#define DIMENSION 2
#endif
#ifdef SUMSIM
#define SYSTEM nKPZ::System
#else
#define SYSTEM nKPZ::SlopeSystem
#endif
#ifdef RLOG
#define DIM nKPZ::LogDim
#else
#define DIM nKPZ::IntDim
#endif

template<class System>
void doMcs(int n, System& system)
{
	unsigned long long int l = n * system.size();
	for(; l > 0; --l)
		system.update();
}

#if 1
template<template<int N, class Dim> class System, class Dim>
void print(std::ostream& o, const System<2, Dim>& system)
{
 	nKPZ::LatticePoint<2, Dim> l;
	unsigned long long int sum = 0;
	for(l.c[0] = 0; l.c[0] < (int)system.dim(); ++l.c[0])
	{
		for(l.c[1] = 0; l.c[1] < (int)system.dim(); ++l.c[1])
		{
			const int t = system.get(l);
			sum += t;
			o << t << "\t";
		}
		o << '\n';
	}
	o << "Sum " << sum << '\n';
}

template<class Dim>
void print(std::ostream& o, const nKPZ::HeightMap<2, Dim>& system)
{
 	nKPZ::LatticePoint<2, Dim> l;
	for(l.c[0] = 0; l.c[0] < (int)system.dim(); ++l.c[0])
	{
		for(l.c[1] = 0; l.c[1] < (int)system.dim(); ++l.c[1])
		{
			o << system.get(l) << "\t";
		}
		o << '\n';
	}
}
#endif

template<template<int N, class Dim> class System, int N, class Dim>
double roughness(const System<N, Dim>& system);

template<int N, class Dim>
double roughness(const nKPZ::System<N, Dim>& system)
{
	nKPZ::HeightMap<N, Dim> h(system.dim());
	h.fromSystem(system);
	//print(std::cout<<"height\n", h);
	return h.roughness();
}

template<int N, class Dim>
double roughness(const nKPZ::SlopeSystem<N, Dim>& system)
{
	return system.roughness();
}

template<template<int N, class Dim> class System, int N, class Dim>
int sim(int argc, char* argv[])
{
	timeAinit;
	{
		timeAstart;
		const int seed = time(0);
		srand(seed);
		std::cout << "Random seed " << seed << '\n';
		dsfmt_gv_init_gen_rand(seed);
		for( int a = 0; a < 10000; ++a) // arbitray large number 
		{
			dsfmt_gv_genrand_close_open();
		}
		timeAstop("initRng");
	}
	if(argc < 2)
		std::cerr << "Usage: scaling mcs [dim]\n";
	int mcs = atoi(argv[1]), dim = 8;
	if(argc > 2)
		dim = atoi(argv[2]);
	double measumentDensity = .1;

	System<N, Dim> system(dim);
	std::cout << "Doing " << mcs << " mcs on a " << N << " dimensional system of size " << (int)system.dim() << std::endl;
	system.printMemStat(std::cout);
	timeAstart
	system.init();
	timeAstop("initSystem");
	timeAstart;
	std::cout << "0\t" << roughness(system) << std::endl;
	timeAstop("calcRoughness");
	//print(std::cout << "system\n", system);

	for(int a = 0; a < mcs;)
	{
		const int nextEval = (int)ceil(exp(log(a+10)+measumentDensity));
		const int nextMcs = std::min(nextEval, mcs) - a;
		a += nextMcs;
		std::cout << a << '\t';
		timeAstart;
		doMcs(nextMcs, system);
		timeAstopS( "mcs_" << nextMcs );
		timeAstart;
		std::cout << roughness(system) << std::endl;
		timeAstop("calcRoughness");
		//print(std::cout << "system\n", system);
	}

}

int main(int argc, char* argv[])
{
	return sim<SYSTEM, DIMENSION, DIM>(argc, argv);
}
