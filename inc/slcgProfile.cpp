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

#include "kmcRandom.h"

#include <iostream>
#include <cmath>
#include <ctime>
#include <cstring>

int main(int argc, char* argv[])
{
	unsigned long long rnd = time(0);
	std::cout << "seed " << rnd << std::endl;
	for(int a = 0; a < 20; ++a)
	{
		rnd = SLCGen(rnd);
		std::cout << (rnd&1) << '\n';
	}
	// bit 0 oscillates with period 2
	// lower bit seem to oscillate, period inceasing with index
	std::cout << "End of bit 0 samples." << std::endl;
	const unsigned int N = 1<<30;
	const int BITS = 32;
	const int JUMPS = 12;
	unsigned int c[BITS];
	unsigned int j[JUMPS], j2[JUMPS];
	memset(c, 0, BITS*sizeof(unsigned int));
	memset(j, 0, JUMPS*sizeof(unsigned int));
	memset(j2, 0, JUMPS*sizeof(unsigned int));
	for(unsigned int a = 0; a < N; ++a)
	{
		rnd = SLCGen(rnd);
		unsigned int t = (rnd ^ (rnd>>32));
		float f = (t >>24)/256.f*JUMPS;
		++j[(int)f];
		f = (t&KMC_M_LCG_RAND_RED)/KMC_LCG_RAND_SUP_RED *JUMPS;
		++j2[(int)f];
		for(int b = 0; t!=0; ++b, t>>=1)
		{
			c[b] += t&1ull;
		}
	}

	std::cout << N << '\t' << sqrt(N) << "\nbits:\n";
	for(int a = 0; a < BITS; ++a)
	{
		std::cout << a << '\t' << c[a] << '\n';
	}
	std::cout << "jumps:\n";
	for(int a = 0; a < JUMPS; ++a)
	{
		std::cout << a << '\t' << j[a] << '\t' << j2[a] << '\n';
	}
	return 0;
}
