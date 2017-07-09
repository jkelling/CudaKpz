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

#include <arghelper.h>

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <cstring>
#include <list>
#include <algorithm>

struct Point
{
	double x, y, dy;

		Point(double x, double y, double dy = NAN) : x(x), y(y), dy(dy) {}

	double effExp(const Point& p) {
		return (log(y) - log(p.y))/(log(x) - log(p.x));
	}
	double effExpErr(const Point& p) {
		return (dy/y + p.dy/p.y)/(log(x) - log(p.x));
	}

	bool operator<(const Point& p) {return x < p.x;}
	bool operator<(const double p) {return x < p;}
};

const int PAIR_TIME = 1;
const int PAIR_IDX = 2;
const int PAIR_IDX_DELTA = 4;
const int PAIR_DENSELY = 8;

template<class Container>
bool advance(Container& c, typename Container::iterator& i, int n)
{
	if(n > 0)
	{
		for(; n > 0 && i != c.end(); --n, ++i);
	}
	else if(n < 0)
	{
		for(; n < 0 && i != c.end(); ++n, --i);
		if(n != 0)
			return false;
	}
	if(n != 0)
		return false;
	return true;
}

int main(int argc, char* argv[])
{
	double minXFrac = 2.;
	double offsetY = 0.;
	int colX = 0, colY = 1, colDY = -1;
	bool verbose = false;
	int pairingMode = PAIR_IDX;
	std::istream &in = std::cin;
	std::ostream &out = std::cout;

	if(minXFrac <= 1)
		return 1;

	for(int i = 1; i < argc; ++i)
	{
		if(!strcmp("-f",argv[i]))
		{
			if(argc <= i+2)
			{
				std::cerr << "Usage: -f colX colY\n";
				return 1;
			}
			colX = atoi(argv[++i]);
			colY = atoi(argv[++i]);
		}
		else if(!strcmp("-e",argv[i]))
		{
			if(argc <= i+1)
			{
				std::cerr << "Usage: -e colDY\n";
				return 1;
			}
			colDY = atoi(argv[++i]);
		}
		else if(!strcmp("-o",argv[i]))
		{
			if(!getArg(offsetY, ++i, argc, argv))
			{
				std::cerr << "Usage: -o offsetY\n";
				return 1;
			}
		}
		else if(!strcmp("-pp",argv[i]))
		{
			if(!getArg(minXFrac, ++i, argc, argv))
			{
				std::cerr << "Usage: -pp paringParam\n";
				return 1;
			}
		}
		else if(!strcmp("-p",argv[i]))
		{
			if(argc <= i+1)
			{
				std::cerr << "Usage: -p pairingMode\n";
				return 1;
			}
			++i;
			if(!strcmp("time",argv[i]))
				pairingMode = PAIR_TIME;
			else if(!strcmp("timeDensely",argv[i]))
				pairingMode = PAIR_TIME | PAIR_DENSELY;
			else if(!strcmp("idx",argv[i]))
				pairingMode = PAIR_IDX;
			else if(!strcmp("idxDelta",argv[i]))
				pairingMode = PAIR_IDX_DELTA;
			else
			{
				std::cerr << "Invalid pairing mode: " << argv[i] << '\n';
				return 1;
			}
		}
		else if(!strcmp("-v",argv[i]))
			verbose = true;
		else return 1;
	}

	int	maxCol = std::max(std::max(colX,colY),colDY);

	std::list<Point> list;
	std::string line;
	while(std::getline(in, line))
	{
		std::istringstream ss(line);
		char t;
		if((!(ss >> t)) || t == '#')
			continue;
		ss.putback(t);
		int c = 0;
		double x, y, dy=NAN, v;
		for(;c <= maxCol; ++c)
		{
			if(!(ss >> v))
			{
				if(ss.bad() || ss.eof())
				{
					c = -1;
					break;
				}
				ss >> line;
				if(c==colX || c==colY || c==colDY)
				{
					std::cerr << "Unreadable line.\n";
					break;
				}
			}
			if(c == colX)
				x = v;
			else if(c == colY)
				y = v-offsetY;
			else if(c == colDY)
				dy = v;
		}
		if(c <= maxCol)
			continue;
		list.push_back(Point(x,y,dy));
	}
	list.sort();
	for(auto a = list.begin(); a != list.end() && a->x == 0; a = list.erase(a));

	auto printEffLine = [&](std::list<Point>::iterator tr, std::list<Point>::iterator ta) -> void {
		out << (ta->x+tr->x)/2. << '\t';
		if(verbose)
			out << ta->x << '\t' << tr->x << '\t';
		out << ta->effExp(*tr);
		if(verbose)
			out << '\t' << ta->y;
		if(colDY >= 0)
			out << '\t' << ta->effExpErr(*tr);
		out << '\n';
	};

	switch (pairingMode)
	{
		case PAIR_TIME | PAIR_DENSELY:
		case PAIR_TIME:
			{
				auto t = list.begin(); ++t;
				for(auto tp = list.begin(); t != list.end(); ++tp)
				{
					while(t->x < tp->x*minXFrac && t != list.end())
						++t;
					if(t == list.end())
						break;
					printEffLine(tp, t);
					if(pairingMode & PAIR_DENSELY)
					{
						auto tt = tp;
						++tt;
						for(++t; t != list.end() && t->x < tt->x*minXFrac; ++t)
							printEffLine(tp, t);
					}
				}
			}
			break;
		case PAIR_IDX:
			{
				int iFrac = ceil(minXFrac);

				auto tr = list.begin();
				auto ta = tr;
				if(!advance(list, ta, iFrac-1))
				{
					std::cerr << "File contains too few points for iFrac= " << iFrac << '\n';
					return 2;
				}

				for(;ta != list.end(); ++tr, advance(list, ta, iFrac))
				{
					printEffLine(tr, ta);
				}
			}
			break;
		case PAIR_IDX_DELTA:
			{
				int iDelta = round(minXFrac);

				auto tr = list.begin();
				auto ta = tr;
				if(!advance(list, ta, iDelta))
				{
					std::cerr << "File contains too few points for iDelta= " << iDelta << '\n';
					return 2;
				}

				for(;ta != list.end(); ++tr, ++ta)
				{
					printEffLine(tr, ta);
				}
			}
			break;
	}

	return 0;
}
