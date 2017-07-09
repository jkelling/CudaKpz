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

#include "filehelper.h"

#include <fstream>
#include <string>
#include <sstream>

int getFiletype(const char* filename)
{
	/*Examine head of file to determine format.*/
	std::ifstream in(filename);
	int lx;
	in >> lx;
	if(in.fail())
		return FT_INVALID;
	in.ignore(10000,'\n');
	std::string s;
	std::getline(in, s);
	if(!s.length())
	{
		//seems to be xyz file
		in.close();
		return FT_XYZ;
	}
	std::istringstream i(s);
	i >> lx;
	if(!in.gcount())
		return FT_INVALID;
	in >> lx;
	if(!in.gcount())
		return FT_INVALID;
	in.close();
	return FT_BIT;
}

void collisionFreeFilename(std::string& orig, const char* suffix)
{
	/* Code harvested from CudaKmc::ConvertXYZ */
	std::ostringstream name;
	std::string localSuffix;
	if(suffix)
		name << orig;
	else
		name << orig.substr(0, orig.find_last_of('.'));
	int offset = name.tellp();
	if(!suffix)
	{
		localSuffix = orig.substr(offset);
		suffix = localSuffix.c_str();
	}
	std::ifstream i;
	name << suffix;
	i.open(name.str().c_str());
	if(i)
	{
		name.seekp(offset);
		name << '_';
		offset = name.tellp();
		int index = 0;
		do {
			i.close();
			++index;
			name.seekp(offset);
			name << index << suffix;
			i.open(name.str().c_str());
		} while(i);
	}
	i.close();
	orig = name.str();
}
