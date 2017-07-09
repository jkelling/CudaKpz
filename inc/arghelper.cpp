/***************************************************************************
*   Copyright (C) 2009, 2013 Jeffrey Kelling <kelling.jeffrey@ages-skripte.org>
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

#include "arghelper.h"

#include <iostream>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <cmath>

bool getArg(std::string& val, int i, int argc, char* argv[])
{
	if(i>= argc)
	{
		std::cerr << "Missing argument.\n";
		return false;
	}
	val = argv[i];
	return true;
}

bool getArg(int& val, int i, int argc, char* argv[])
{
	if(i>= argc)
	{
		std::cerr << "Missing argument.\n";
		return false;
	}
	std::istringstream ss(argv[i]);
	int a;
	ss >> a;
	if(!ss.eof())
	{
		std::cerr << "Expected an integral number.\n";
		return false;
	}
	val = a;
	return true;
}

bool getArg(long long int& val, int i, int argc, char* argv[])
{
	if(i>= argc)
	{
		std::cerr << "Missing argument.\n";
		return false;
	}
	std::istringstream ss(argv[i]);
	long long int a;
	ss >> a;
	if(!ss.eof())
	{
		std::cerr << "Expected an integral number.\n";
		return false;
	}
	val = a;
	return true;
}

bool getArg(unsigned int& val, int i, int argc, char* argv[])
{
	if(i>= argc)
	{
		std::cerr << "Missing argument.\n";
		return false;
	}
	std::istringstream ss(argv[i]);
	int a;
	ss >> a;
	if(!ss.eof())
	{
		std::cerr << "Expected an unsigned integral number.\n";
		return false;
	}
	val = a;
	return true;
}

bool getArg(long double& val, int i, int argc, char* argv[], int flags)
{
	if(i>= argc)
	{
		std::cerr << "Missing argument.\n";
		return false;
	}
	std::istringstream ss(argv[i]);
	long double a;
	ss >> a;
	if(!ss.eof())
	{
		if(flags & (getArgFlags::FLOAT_ALLOW_NAN | getArgFlags::FLOAT_ALLOW_INF))
		{
			std::string s(argv[i]);
			std::for_each(s.begin(), s.end(), tolower);
			if((flags & getArgFlags::FLOAT_ALLOW_NAN) && s == "nan")
			{
				val = NAN;
				return true;
			}
			else if((flags & getArgFlags::FLOAT_ALLOW_INF)
					&& s.substr(s.length()-3, s.length()-1) == "inf")
			{
				if(s.length() == 3 || (s.length()==4 && s[0] == '+'))
				{
					val = +INFINITY;
					return true;
				}
				else if(s.length()==4 && s[0] == '-')
				{
					val = -INFINITY;
					return true;
				}
			}
		}
		std::cerr << "Expected a floatingpoint number.\n";
		return false;
	}
	val = a;
	return true;
}

bool getArg(double& val, int i, int argc, char* argv[], int flags)
{
	if(i>= argc)
	{
		std::cerr << "Missing argument.\n";
		return false;
	}
	std::istringstream ss(argv[i]);
	double a;
	ss >> a;
	if(!ss.eof())
	{
		if(flags & (getArgFlags::FLOAT_ALLOW_NAN | getArgFlags::FLOAT_ALLOW_INF))
		{
			std::string s(argv[i]);
			std::for_each(s.begin(), s.end(), tolower);
			if((flags & getArgFlags::FLOAT_ALLOW_NAN) && s == "nan")
			{
				val = NAN;
				return true;
			}
			else if((flags & getArgFlags::FLOAT_ALLOW_INF)
					&& s.substr(s.length()-3, s.length()-1) == "inf")
			{
				if(s.length() == 3 || (s.length()==4 && s[0] == '+'))
				{
					val = +INFINITY;
					return true;
				}
				else if(s.length()==4 && s[0] == '-')
				{
					val = -INFINITY;
					return true;
				}
			}
		}
		std::cerr << "Expected a floatingpoint number.\n";
		return false;
	}
	val = a;
	return true;
}

bool getArg(float& val, int i, int argc, char* argv[], int flags)
{
	if(i>= argc)
	{
		std::cerr << "Missing argument.\n";
		return false;
	}
	std::istringstream ss(argv[i]);
	long double a;
	ss >> a;
	if(!ss.eof())
	{
		if(flags & (getArgFlags::FLOAT_ALLOW_NAN | getArgFlags::FLOAT_ALLOW_INF))
		{
			std::string s(argv[i]);
			std::for_each(s.begin(), s.end(), tolower);
			if((flags & getArgFlags::FLOAT_ALLOW_NAN) && s == "nan")
			{
				val = NAN;
				return true;
			}
			else if((flags & getArgFlags::FLOAT_ALLOW_INF)
					&& s.substr(s.length()-3, s.length()-1) == "inf")
			{
				if(s.length() == 3 || (s.length()==4 && s[0] == '+'))
				{
					val = +INFINITY;
					return true;
				}
				else if(s.length()==4 && s[0] == '-')
				{
					val = -INFINITY;
					return true;
				}
			}
		}
		std::cerr << "Expected a floatingpoint number.\n";
		return false;
	}
	val = a;
	return true;
}

bool getArg(const char*& val, int i, int argc, char* argv[])
{
	if(i>= argc)
	{
		std::cerr << "Missing argument.\n";
		return false;
	}
	val = argv[i];
	return true;
}
