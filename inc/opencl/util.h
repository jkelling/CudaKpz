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

#ifndef KMC_OPENCL_UTIL_H
#define KMC_OPENCL_UTIL_H

#include <CL/cl.h>
#include <string>
#include <iostream>
#include <sstream>

namespace KmcCl
{
#define checkErr(err, name) \
{ \
	if (err != CL_SUCCESS) \
	{ \
		std::cerr << "ERROR: " << name << " (" << err << ')' << " in " << __FILE__ << " at line " << __LINE__ << std::endl; \
		exit(EXIT_FAILURE); \
	} \
}

class CuClParser
{
	private:
	void globalConstExpr(std::ostream& o, const std::string& line, std::istream& file);
	bool subLoadPatchCode(std::string dir, const std::string& filename, std::ostream& o);
		
	protected:

	std::ostringstream m_code;

	virtual bool paste(std::ostream& o, const std::string& directive) {return false;}
	public:
		CuClParser() {}
	virtual ~CuClParser() {}

	int loadPatchCode(const std::string &file);
	cl_program makeProgram(cl_context context, cl_device_id* device);

	const std::ostringstream& code() const {return m_code;}

	protected:

	struct ScopeState
	{
		int scope;
		int ignoreState;

		static const int NONE = 0;
		static const int STRING = 1;
		static const int C_COMMENT = 2;
		static const int MACRO = 3;

		void parse(const std::string& line);

			ScopeState () : scope(0), ignoreState(0) {}
		inline bool global() const {return scope == 0 && ignoreState == NONE;}
	};
};

void buildProgram(cl_program& program, cl_device_id* device);
}

#endif
