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
#include <sstream>
#include <list>
#include <vector>
#include <string>
#include <cmath>
#include <cstring>

class Data
{
	typedef std::vector<double> Row;

	std::list<Row> m_data;
	int m_longestRow;
	int m_colX;
	double m_min, m_max;

	friend std::ostream& operator<<(std::ostream& o, const Row& d);

	void prebin();
	void bin(const double width);

	public:

		Data(int colX = 0) : m_longestRow(2), m_colX(colX) {}

	void read(std::istream& i);
	int readLine(std::istream& i);

	double binFixNumber(int n);
	int binFixWidth(double w);

	friend std::ostream& operator<<(std::ostream& o, const Data& d);
};


int main(int argc, char* argv[])
{
	double width = NAN;
	int bins = 8;

	if(argc > 2)
	{
		if(!strcmp("-w", argv[1]))
		{
			width = atof(argv[2]);
			if(width <= 0)
			{
				std::cerr << "Expected width > 0\n";
				return 1;
			}
		}
		else if(!strcmp("-n", argv[1]))
		{
			bins = atoi(argv[2]);
			if(bins < 1)
			{
				std::cerr << "Expected bins >= 1\n";
				return 1;
			}
		}
		else
		{
			std::cerr << "Invalid argument: " << argv[1] << '\n';
			return 1;
		}
	}
	else if(argc > 1)
	{
		std::cerr << "Invalid argument: " << argv[1] << " (possibly missing parameter)\n";
		return 1;
	}

	Data data;
	data.read(std::cin);
	if(std::isnan(width))
	{
		const double w = data.binFixNumber(bins);
		std::cout << "#bin width " << w << '\n';
	}
	else
	{
		const int n = data.binFixWidth(width);
		std::cout << "#bins " << n << '\n';
	}
	std::cout << data;

	return 0;
}

int Data::binFixWidth(double w)
{
	if(w <= 0)
		return 0;
	prebin();
	bin(w);
	return (m_max-m_min)/w;
}

double Data::binFixNumber(int n)
{
	if(n < 1)
		return NAN;
	prebin();
	const double bw = (m_max-m_min)/n;
	bin(bw);
	return bw;
}

void Data::prebin()
{
	for(std::list<Row>::iterator a = m_data.begin(); a != m_data.end(); ++a)
		if(a->size() <= m_colX)
		{
			a = m_data.erase(a);
			if(a == m_data.end())
				break;
		}

	m_data.sort([this](const Row& a, const Row& b){return a[m_colX] < b[m_colX];});

	m_min = m_data.front()[m_colX];
	m_max = m_data.back()[m_colX];
}

void Data::bin(const double width)
{
	std::list<Row> binned;
	for(double m = m_min+width; m <= m_max; m += width)
	{
		Row row;
		row.resize(m_longestRow);
		row[m_colX] = m-width/2.;
		for(int s = 0; s < m_colX; ++s)
			row[s] = 0;
		for(int s = m_colX+1; s < row.size(); ++s)
			row[s] = 0;
		for(std::list<Row>::iterator a = m_data.begin(); a != m_data.end() && (*a)[m_colX] < m; a = m_data.erase(a))
		{
			for(int s = 0; s < m_colX; ++s)
			{
				const double t = (*a)[s];
				if(!std::isnan(t))
					row[s] += t;
			}
			for(int s = m_colX+1; s < a->size(); ++s)
			{
				const double t = (*a)[s];
				if(!std::isnan(t))
					row[s] += t;
			}
		}
		binned.push_back(row);
	}
	m_data.swap(binned);
}

void Data::read(std::istream& i)
{
	std::string s;
	while(getline(i,s))
	{
		std::istringstream il(s);
		readLine(il);
	}
}

int Data::readLine(std::istream& i)
{
	{
		char t;
		i >> t;
		if(t == '#')
			return 0;
		else
			i.putback(t);
	}
	Row row;
	row.resize(m_longestRow);
	double t;
	int a = 0;
	while(i >> t)
	{
		if(a >= m_longestRow)
		{
			row.resize(a+1);
			m_longestRow = a+1;
		}
		row[a] = t;
		++a;
	}
	if(a)
		m_data.push_back(row);
	return a;
}

std::ostream& operator<<(std::ostream& o, const Data::Row& d)
{
	if(d.size())
	{
		o << d[0];
		for(int a = 1; a < d.size(); ++a)
		{
			o << '\t' << d[a];
		}
	}
	return o;
}

std::ostream& operator<<(std::ostream& o, const Data& d)
{
	for(std::list<Data::Row>::const_iterator a = d.m_data.begin(); a != d.m_data.end(); ++a)
	{
		o << *a << '\n';
	}
	return o;
}
