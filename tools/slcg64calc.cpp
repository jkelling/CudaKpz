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

#include <iostream>
#include <cstdlib>

#define KMC_L_LCG_N 4
#include <kmcRandom.h>

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		std::cerr << "Usage slcg64calc num_generators [comment [...]]\n";
		return 0;
	}
	int a = atoi(argv[1]);
	if(a <= 0)
	{
		std::cerr << "invalid input: " << a << '\n';
		return 1;
	}
	unsigned long long askip = SLCGskipA(a);
	unsigned long long cskip = SLCGskipC(a);

	if(argc > 2)
	{
		std::cout << "\n// Instatiation for";
		for(int i = 2; i < argc; ++i)
			std::cout << ' ' << argv[i];
	}
	std::cout << "\ntemplate<>\nstruct SLCGskipAt<" << a
		<< "> {\n\tstatic const unsigned long long a = "
		<< askip << "ull;\n};\n\ntemplate<>\nstruct SLCGskipCt<" << a
		<< "> {\n\tstatic const unsigned long long c = "
		<< cskip << "ull;\n};\n";

	return 0;
}
