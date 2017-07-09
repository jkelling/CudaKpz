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

#ifndef KMC_BENCHMARKING_H
#define KMC_BENCHMARKING_H

#include <sys/time.h>
#include <sstream>
#include <iostream>

// always-on
#define timeInit_P(marker) \
		timeval t ## marker ## 1, t ## marker ## 2;
#define timeStart_P(marker) \
		gettimeofday(&t ## marker ## 1,0); 
#define timeStop_P(marker,arg) \
		gettimeofday(&t ## marker ## 2,0); \
		{ \
			std::ostringstream s; \
			s << "benchmark-" arg "\t" << std::fixed << \
				((float)(t ## marker ## 2 . tv_sec - t ## marker ## 1 . tv_sec) \
				+ ((float)(t ## marker ## 2 . tv_usec - t ## marker ## 1 . tv_usec))*1e-6) \
				<< '\n'; \
			std::cerr << s.str(); \
		}
#define timeStopS_P(marker,arg) \
		gettimeofday(&t ## marker ## 2,0); \
		{ \
			std::ostringstream s; \
			s << "benchmark-" << arg << '\t' << std::fixed << \
				((float)(t ## marker ## 2 . tv_sec - t ## marker ## 1 . tv_sec) \
				+ ((float)(t ## marker ## 2 . tv_usec - t ## marker ## 1 . tv_usec))*1e-6) \
				<< '\n'; \
			std::cerr << s.str(); \
		}

// switchable
#ifdef KMC_SWITCH_BENCHMARK
#define timeAinit \
		timeval tA1,tA2;
#define timeAstart \
		gettimeofday(&tA1,0); 
#define timeAstop(arg) \
		gettimeofday(&tA2,0); \
		{ \
			std::ostringstream s; \
			s << "benchmark-" arg "\t" << std::fixed << \
				((float)(tA2 . tv_sec - tA1 . tv_sec) \
				+ ((float)(tA2 . tv_usec - tA1 . tv_usec))*1e-6) \
				<< '\n'; \
			std::cerr << s.str(); \
		}
#define timeAstopS(arg) \
		gettimeofday(&tA2,0); \
		{ \
			std::ostringstream s; \
			s << "benchmark-" << arg << '\t' << std::fixed << \
				((float)(tA2 . tv_sec - tA1 . tv_sec) \
				+ ((float)(tA2 . tv_usec - tA1 . tv_usec))*1e-6) \
				<< '\n'; \
			std::cerr << s.str(); \
		}

#define timeInit(marker) \
		timeval t ## marker ## 1, t ## marker ## 2;
#define timeStart(marker) \
		gettimeofday(&t ## marker ## 1,0); 
#define timeStop(marker,arg) \
		gettimeofday(&t ## marker ## 2,0); \
		{ \
			std::ostringstream s; \
			s << "benchmark-" arg "\t" << std::fixed << \
				((float)(t ## marker ## 2 . tv_sec - t ## marker ## 1 . tv_sec) \
				+ ((float)(t ## marker ## 2 . tv_usec - t ## marker ## 1 . tv_usec))*1e-6) \
				<< '\n'; \
			std::cerr << s.str(); \
		}
#define timeStopS(marker,arg) \
		gettimeofday(&t ## marker ## 2,0); \
		{ \
			std::ostringstream s; \
			s << "benchmark-" << arg << '\t' << std::fixed << \
				((float)(t ## marker ## 2 . tv_sec - t ## marker ## 1 . tv_sec) \
				+ ((float)(t ## marker ## 2 . tv_usec - t ## marker ## 1 . tv_usec))*1e-6) \
				<< '\n'; \
			std::cerr << s.str(); \
		}
#else
#define timeAinit
#define timeAstart
#define timeAstop(str)
#define timeAstopS(str)
#define timeInit(marker)
#define timeStart(marker)
#define timeStop(marker,str)
#define timeStopS(marker,str)
#endif

#endif
