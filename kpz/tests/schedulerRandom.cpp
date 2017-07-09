/***************************************************************************
*   Copyright 2013 Jeffrey Kelling <j.kelling@hzdr.de>
*                  Helmholtz-Zentrum Dresden-Rossendorf
*                  Institute of Ion Beam Physics and Materials Research
*
*   This program is free software; you can redistribute it and/or
*   modify it under the terms of the GNU General Public
*   License as published by the Free Software Foundation; either
*   version 2 of the License, or (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program; if not, write to the Free Software
*   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
***************************************************************************/

#include "../systemSize.h"
#include "../schedulerService.h"

#include <iostream>

int main()
{
	Kpz::SchedulerService<Kpz::SystemSizeCPU> schedulerA(5,5), schedulerB(5,5);

	static const int N = 10;

	for (int a = 0; a < N; ++a)
		std::cout << "rand\t" << a << ":\t" << schedulerA.genrand_close_open()
			<< '\t' << schedulerB.genrand_close_open() << '\n';

	std::cout << "-- sync --\n";
	schedulerB.copyDsfmt(schedulerA.dsfmt());

	for (int a = 0; a < N; ++a)
		std::cout << "rand\t" << a << ":\t" << schedulerA.genrand_close_open()
			<< '\t' << schedulerB.genrand_close_open() << '\n';

	return 0;
}
