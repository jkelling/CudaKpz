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

#include <cstring>
#include <cmath>

#include <arghelper.h>
#include <combineTables.h>

using namespace Kmc;

int main (int argc, char* argv[])
{
	int idxX = 0;
	std::list<int> idxY{1};
	double value = NAN;
	bool printFilename = true;
	std::ostream& out(std::cout);
	
	for(int i = 1; i < argc; ++i)
	{
		if(!strcmp("-x", argv[i]))
		{
			if(!getArg(idxX, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("-y", argv[i]))
		{
			std::string tmp;
			if(!getArg(tmp, ++i, argc, argv))
				return 1;
			if(!Kmc::parseIdxYList(idxY, tmp))
				return 1;
		}
		else if(!strcmp("-v", argv[i]) || !strcmp("--verbose", argv[i]))
		{
			Kmc::WARN.on();
		}
		else if(!strcmp("--noPrintFilename", argv[i]))
			printFilename = false;
		else if(!strcmp("-s", argv[i]))
		{
			if(!getArg(value, ++i, argc, argv))
				return 1;
		}
		else if(std::isnan(value))
		{
			if(!getArg(value, i, argc, argv))
				return 1;
		}
		else
		{
			Kmc::TableFileView file(idxX, idxY, argv[i]);
			if(!file.open())
			{
				std::cerr << "[EE] Failed to open file " << file.name() << '\n';
				continue;
			}
			Kmc::TableFileView::Pair prev = file.currentPair();
			// assuming ascending ordering
			while(file.good())
			{
				if(value < file.x()) // found
				{
					if(printFilename)
						out << file.name() << '\t';

					const double c = (value-prev.x()) / (file.x() - prev.x());
					for(int a = 0; a < file.nY(); ++a)
					{
						// linear interpolation
						const double y = prev.y(a) + (file.y(a) - prev.y(a))*c;
						out << ((a==0)? "" : "\t") << y;
					}
					out << '\n';
					break;
				}
				prev = file.currentPair();
				file.nextLine();
			}
			if(file.good())
				WARN() << "Value not found in file " << file.name() << '\n';
		}
	}

	return 0;
}
