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


template<class Dim, class Elem, class Output>
void Kmc::writePM3D(const Map<2, Dim, Elem> &map, std::ostream& o, double cellDim, Output out)
{
	size_t index = 0;
	for(int y = 0; y < (int)map.dim(1); ++y)
	{
		for(int x = 0; x < (int)map.dim(0); ++x)
		{
			o << x*cellDim << '\t' << y*cellDim << '\t';
			out(o, map[index]) << '\n';
			++index;
		}
		o << '\n';
	}
}

template<int N, class Dim, class Elem, class Output = std::ostream& (std::ostream&, const Elem&)>
void Kmc::writeASCII(const Map<N, Dim, Elem> &map, std::ostream& o, double cellDim = 1., Output out = outOperator<Elem>)
{
	typename Map<2, Dim, Elem>::Point p(0);
	p.coord[0] = -1;
	for(int a = 0; a < map.size(); ++a)
	{
		int s;
		for(s = 0; s < N; ++s)
		{
			++p.coord[s];
			if(p.coord[s] >= (int)map.dim(s))
			{
				p.coord[s] = 0;
				o << cellDim*p.coord[s] << '\t';
			}
			else
				break;
		}
		for(; s < N; ++s)
			o << cellDim*p.coord[s] << '\t';
		out(o, map[a]) << '\n';
	}
}

template<class Dim, class Elem, class Output>
void Kmc::writeOctaveMat(const Map<2, Dim, Elem> &map, std::ostream& o, const char* name, Output out)
{
	o << "# Kmc::Map<2>\n# name: " << name << "\n# type: matrix\n# rows: " << (int)map.dim(1)
		<< "\n# columns: " << (int)map.dim(0) << '\n';
	size_t index = 0;
	for(int y = 0; y < (int)map.dim(1); ++y)
	{
		for(int x = 0; x < (int)map.dim(0); ++x)
		{
			out(o << ' ', map[index]);
			++index;
		}
		o << '\n';
	}
}

template<int N, class Dim, class Elem>
Elem Kmc::variance(const Kmc::Map<N, Dim, Elem> &map)
{
	Elem sum = 0, sum2 = 0;
	for(size_t a = 0; a < map.size(); ++a)
	{
		sum += map[a];
		sum2 += map[a]*map[a];
	}
	return (sum2-sum*sum/map.size())/(map.size()-1);
}
