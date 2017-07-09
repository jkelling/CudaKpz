/***************************************************************************
*   Copyright 2010 - 2012 Jeffrey Kelling <j.kelling@hzdr.de>
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

/*
 * These functions are only compatible with POSIX, not with Windows.
 */

#include <string>

/*! \ingroup DGio
 * \{
 * \param s filename
 * \return basename of path given in \p s
 */
std::string basename(const std::string& s)
{
	std::string::size_type begin = s.find_last_of('/');
	//if(begin == std::string::npos)
	//	begin = -1; //npos == -1
	++begin;
	//basename is what comes after the last '/', ore everything if none present
	return s.substr(begin);
}

/*! \param s filename
 * \return dirname of path given in \p s
 */
std::string dirname(const std::string& s)
{
	std::string::size_type end = s.find_last_of('/');
	if(end == std::string::npos)
		return "";
	//dirname is everything which is not basename
	return s.substr(0,++end);
}

#include <iostream>
/*!
 * Changes \p d to the path the location given by \p s relative to \p d.
 * \warning Assumes \p d to have a trailing '/'.
 * Treats \p s as dirname(\p s).
 * \param d basepath (will be modified)
 * \param s relative path
 * \return reference to d
 */
std::string& updateDirname(std::string& d, const std::string& s)
{
	if(!s.size())
		return d;
	if((s[0]=='/') || !d.size())
		return d = s.substr(0,s.find_last_of('/')+1);
	std::string::size_type pos = 0, newpos;
	while((newpos = s.find_first_of('/', pos)) != std::string::npos)
	{
		std::string tmp(s.substr(pos, ++newpos-pos));
		if(tmp == "./")
			continue;
		if(tmp == "../")
		{
			if(!d.size())
			{
				d="../";
			}
			else if(d.size()==1)
			{
				if(d[0] == '/')
					return d;
				d="";
			}
			else
			{	
				d.erase(d.size()-1);
				std::string::size_type index = d.find_last_of('/');

				if(index == std::string::npos)
				{
					if(d == "..")
						d += "/../";
					else
						d = "";
				}
				else if(d.substr(index) == "/..")
					d += "/../";
				else
				{
					//if we are dealing with the last remainder of a relative path dir is empty afterwards, since npos+1==0
					d.erase(index+1);
				}
			}
		}
		else
		{
			d += tmp;
		}
		pos = newpos;
	}
	return d;
}
