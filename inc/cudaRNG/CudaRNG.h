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

#ifndef KMC_CUDA_MT_H
#define KMC_CUDA_MT_H

#include "MersenneTwister_NVGPUSDK.h"

class CudaRNG
{
	unsigned int* m_dRandom;
	int m_nPerRng;

	public:

		CudaRNG();
		CudaRNG(const char* file);
		~CudaRNG();

	inline unsigned int* dRandom() {return m_dRandom;}
	inline int n() {return m_nPerRng*RNG_COUNT;}

	void setCount (int n);

	void generate();
	void generateLCG();
	void generateDSFMT();
	void generateLCGS();
	void generateRC4();
	void generateTwofish();
	bool twofishTest();

	static const int THREADS = 512; // Twofish, LCG
	static const int THREADS_RED = 128; // MT, LCGS
	static const int L_MT_N_PER_RNG_MUL = 2; // MT, LCGS
	static const int L_RC4_N_PER_RNG_MUL = 3; // RC4
	static const int BLOCKS = 32;
	static const int RNG_COUNT = THREADS*BLOCKS; // max eg. Twofish, LCG

};

#endif
