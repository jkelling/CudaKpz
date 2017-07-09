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

#include "util.h"

#include <fstream>
#include <cstring>
#include <unistd.h>

#include <relativePath.h>

void KmcCl::CuClParser::globalConstExpr(std::ostream& o, const std::string& line, std::istream& file)
{
	/*
	This function implements a quick-and-dirty way to replace constant declarations by equivalent #defines
	*/
	o << "#define ";
	std::string last, tmp;
	std::istringstream i(line);
	while(i>>tmp)
	{
		if(tmp[0] == '=')
		{
			o << last << " (";
			break;
		}
		else
			last.swap(tmp);
	}
	if(tmp.size() > 1)
		last = tmp.substr(2);
	else
		last = "";
	last.append(std::istreambuf_iterator<char>(i), (std::istreambuf_iterator<char>()));
	std::string::size_type pos = last.find_first_of(';');
	o << last.substr(0, pos );
	while(pos == std::string::npos)
	{
		getline(file, last);
		pos = last.find_first_of(';');
		o << last.substr(0, pos );
	}
	o << ")\n";
}

bool KmcCl::CuClParser::subLoadPatchCode(std::string dir, const std::string& filename, std::ostream& o) 
{
	dir = updateDirname(dir, filename);

	std::ifstream f((dir+basename(filename)).c_str());
	if(!f.is_open())
	{
		std::cerr << (filename + " -- does not exist or is not a regular file.").c_str() << '\n';
		return false;
	}
	
	std::string line;
	ScopeState scope;
	int guard = 0, guardLevel = 0;
	static const int GUARD_NONE = 0;
	static const int GUARD_USE = 1;
	static const int GUARD_IGN = 2;
	static const int GUARD_FINAL = 4;

	o << "//BEGIN " << filename << '\n';

	std::string tmp;
	while(getline(f, tmp))
	{
		std::istringstream i(tmp);
		i >> tmp;

		if(guard & GUARD_IGN)
		{
			if(tmp == "#endif")
			{
				--guardLevel;
				if(!guardLevel)
					guard = GUARD_NONE;
			}
			else if(tmp == "#else")
			{
				if(guardLevel == 1 && !(guard & GUARD_FINAL))
					guard = GUARD_USE;
			}
			else if(tmp == "#elif")
			{
				if(guardLevel == 1 && !(guard & GUARD_FINAL))
				{ // connot parse this convert to #if
					guardLevel = 0;
					guard = GUARD_NONE;
					o << "#if " << i.str() << '\n';
				}
			}
			else if(tmp == "#if" || tmp == "#ifdef" || tmp == "#ifndef")
				++guardLevel;
		}
		else if(tmp == "#include") //TODO: protect from mutual inclusion
		{
			char bracket;
			if(!(i>>bracket) || (bracket != '<' && bracket != '"'))
			{
				std::cerr << "Invalid include: "<< bracket << i.str() << " in file " << filename  << '\n';
				return false;
			}

			tmp = i.str().substr(i.tellg());
			if(bracket == '<')//TODO: implement search paths
			{
				std::string::size_type pos = tmp.find_first_of('>');
				tmp = std::string(getenv("KMCSVNDIR")) + "/inc/" + tmp.substr(0, pos);
			}
			else
			{
				std::string::size_type pos = tmp.find_first_of(bracket);
				tmp = tmp.substr(0, pos);
			}
			if(!subLoadPatchCode(dir, tmp, o))
			{
				std::cerr << "... in " << filename << '\n';
				return false;
			}
		}
		else if(tmp == "#pragma")
		{
			if(!(i>>tmp))
			{
				std::cerr << "Empty #pragma\n";
				return false;
			}
			if(tmp == "KMC_OCL_RUNTIME_MOD")
			{
				if(!(i >> tmp))
				{
					std::cerr << "Empty KMC_OCL_RUNTIME_MOD\n";
					return false;
				}
				if(tmp == "paste")
				{
					if(!paste(o, i.str().substr(i.tellg())))
					{
						std::cerr << "Invalid KMC_OCL_RUNTIME_MOD paste" << i.str().substr(i.tellg()) << '\n';
					}
				}
				else
				{
					std::cerr << "Invalid KMC_OCL_RUNTIME_MOD " << tmp << '\n';
				}
			}
			else if(tmp == "unroll")
			{ // remove nvcc #pragma unroll
			}
			else // uninteresting #pragma, put back
			{
				o << i.str() << '\n';
			}
		}
		else if(scope.global() && tmp == "const")
		{
			globalConstExpr(o, i.str().substr(i.tellg()), f);
		}
		else
		{
			if(tmp == "#ifdef")
			{
				if(guard == GUARD_USE)
					++guardLevel;
				else
				{
					i >> tmp;
					if(tmp == "__OPENCL_VERSION__")
					{
						guard = GUARD_USE;
						guardLevel = 1;
						continue;
					}
				}
			}
			else if(tmp == "#ifndef")
			{
				if(guard == GUARD_USE)
					++guardLevel;
				else
				{
					i >> tmp;
					if(tmp == "__OPENCL_VERSION__")
					{
						guard = GUARD_IGN;
						guardLevel = 1;
						continue;
					}
				}
			}
			else if(guard == GUARD_USE)
			{
				if(tmp == "#if")
					++guardLevel;
				else if(tmp == "#endif")
				{
					--guardLevel;
					if(!guardLevel)
					{
						guard = GUARD_NONE;
						continue;
					}
				}
				else if(tmp == "#elif")
				{
					if(guardLevel == 1)
					{
						guard = GUARD_IGN | GUARD_FINAL;
						continue;
					}
				}
				else if(tmp == "#else")
				{
					if(guardLevel == 1)
					{
						guard = GUARD_IGN | GUARD_FINAL;
						continue;
					}
				}
			}
			scope.parse(i.str());
			o << i.str() 
				//<< " //" << scope.scope << '-' << scope.ignoreState
				<< '\n';
		}
	}
	o << "//END " << filename << '\n';
	return true;
}

int KmcCl::CuClParser::loadPatchCode(const std::string &file)
{
	m_code.str("");
	if(!subLoadPatchCode("", file, m_code))
	{
		std::cerr << "Parsingerror.\n";
		exit(1);
	}
	return m_code.str().size();
}

cl_program KmcCl::CuClParser::makeProgram(cl_context context, cl_device_id* device)
{
	const int size = m_code.str().size();
	char* files = new char[size+1];
	memcpy(files, m_code.str().c_str(), size+1);
	cl_int err;
	cl_program prog = clCreateProgramWithSource(context, 1, (const char**)&files, 0, &err);
	checkErr(err,"makeProgramm");
	buildProgram(prog, device);
	delete[] files;
	return prog;
}

void KmcCl::buildProgram(cl_program& program, cl_device_id* device)
{
	cl_int err;
	err = clBuildProgram(program, 1/*devicenumber*/, device/*devices*/, NULL, NULL, NULL);

	if(err != CL_SUCCESS)
	{
		if(err == CL_BUILD_PROGRAM_FAILURE)
		{
			cl_int logStatus;
			char * buildLog = NULL;
			size_t buildLogSize = 0;
			logStatus = clGetProgramBuildInfo(program,
				*device,
				CL_PROGRAM_BUILD_LOG,
				buildLogSize,
				buildLog,
				&buildLogSize);
			checkErr(logStatus, "clGetProgramBuildInfo() failed");

			buildLog = (char*)malloc(buildLogSize);
			if(buildLog == NULL)
			{
				checkErr(err, "Failed to allocate host memory.(buildLog)");
			}

			memset(buildLog, 0, buildLogSize);

			logStatus = clGetProgramBuildInfo(program,
				*device,
				CL_PROGRAM_BUILD_LOG,
				buildLogSize,
				buildLog,
				NULL);
			if(logStatus != CL_SUCCESS)
			{
				free(buildLog);
				checkErr(err, "clGetProgramBuildInfo() failed");
			}

			std::cout << " \n\t\t\tBUILD LOG\n";
			std::cout << " ************************************************\n";
			std::cout << buildLog << std::endl;
			std::cout << " ************************************************\n";
			free(buildLog);
		}
		checkErr(err, "clBuildProgram()");
	}
}

void KmcCl::CuClParser::ScopeState::parse(const std::string& line)
{
	if(ignoreState != MACRO)
	{
		std::istringstream i(line);
		char a;
		i >> a;
		if(a == '#')
			ignoreState = MACRO;
	}
	if(ignoreState == MACRO)
	{
		if(line[line.size()-1] != '\\')
			ignoreState = NONE;
		return;
	}
	for(std::string::const_iterator i = line.begin(); i != line.end(); ++i)
	{
		switch(ignoreState)
		{
		case STRING:
			switch (*i)
			{
				case '"':
					ignoreState = NONE;
					break;
				case '\\':
					++i; // escape charater: ignore next one
					if(i == line.end())
						return;
					break;
			}
			break;
		case C_COMMENT:
			if(*i == '*')
			{
				++i;
				if(i == line.end())
					return;
				if(*i == '/')
					ignoreState = NONE;
			}
			break;
		default:
			switch(*i)
			{
				case '"':
					ignoreState = STRING;
					break;
				case '/':
					++i;
					if(i == line.end())
						return;
					switch(*i)
					{
						case '*':
							ignoreState = C_COMMENT;
							break;
						case '/': // C++-style comment
							if(line[line.size()-1] == '\\')
								ignoreState = MACRO;
							return;
					}
					break;
				case '{':
					++scope;
					break;
				case '}':
					--scope;
			}
		}
	}
}
