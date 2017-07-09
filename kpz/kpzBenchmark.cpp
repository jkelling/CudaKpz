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


#include "scheduler.h"

#include <iostream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <sys/time.h>
#include <iomanip>

int main(int argc, char* argv[])
{
	int x, y;
	if(argc < 2)
		x = y = 10;
	else if(argc < 3)
		x = y = atoi(argv[1]);
	else
	{
		x = atoi(argv[1]);
		y = atoi(argv[2]);
	}
	const int mcs = 100;
	const int rounds = 30;
	Kpz::Scheduler* scheduler = new Kpz::Scheduler(x,y);

	timeval t1, t2;
	double sum, sqsum;
	
	scheduler->cout() << "#mcs=" << mcs
		<< "\n#time\tstdev" << std::endl;

	sum = sqsum = 0;
	for(int r = 0; r < rounds; ++r)
	{
		scheduler->initHomogeneous();
		scheduler->pushSystem();
		scheduler->randomize();
		gettimeofday(&t1, 0);
		scheduler->mcs(mcs);
		gettimeofday(&t2, 0);
		const double val = (double)(t2.tv_sec - t1.tv_sec) + ((double)(t2.tv_usec - t1.tv_usec))*1e-6;
		sum += val;
		sqsum += val*val;
	}
	sum /= rounds;
	sqsum /= rounds;
	scheduler->cout() << sum << '\t' << sqrt(sqsum - sum*sum) << std::endl;

	delete scheduler;
	return 0;
}
