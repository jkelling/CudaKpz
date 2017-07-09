/***************************************************************************
*   Copyright (C) 2009, 2012 - 2013 Jeffrey Kelling <kelling.jeffrey@ages-skripte.org>
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

#ifndef ARGHELPER_H 
#define ARGHELPER_H

#include <string>

/*! \ingroup DGio
 * \{ */

namespace getArgFlags
{
	const int DEFAULT = 0;
	const int FLOAT_ALLOW_NAN = 1;
	const int FLOAT_ALLOW_INF = 2;
}

/*! Parse next argument on the command-line as string.
 * \param val reference to where to write the parsed value
 * \param i index of the next argument
 * \param argc number of arguments
 * \param argv list of arguments
 * \return whether argument was successfully parsed
 * \ingroup DGio
 */
bool getArg(std::string& val, int i, int argc, char* argv[]);
/*! \overload getArg(std::string&,int,int,char*[]) */
bool getArg(int& val, int i, int argc, char* argv[]);
/*! \overload getArg(std::string&,long long int,int,char*[]) */
bool getArg(long long int& val, int i, int argc, char* argv[]);
/*! \overload getArg(std::string&,int,int,char*[]) */
bool getArg(unsigned int& val, int i, int argc, char* argv[]);
/*! \overload getArg(std::string&,int,int,char*[]) */
bool getArg(long double& val, int i, int argc, char* argv[], int flags = getArgFlags::DEFAULT);
/*! \overload getArg(std::string&,int,int,char*[]) */
bool getArg(double& val, int i, int argc, char* argv[], int flags = getArgFlags::DEFAULT);
/*! \overload getArg(std::string&,int,int,char*[]) */
bool getArg(float& val, int i, int argc, char* argv[], int flags = getArgFlags::DEFAULT);
/*! \overload getArg(const char*&,int,int,char*[]) */
bool getArg(const char*& val, int i, int argc, char* argv[]);

//!\}
#endif
