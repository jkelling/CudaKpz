/***************************************************************************
*   Copyright 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef KMC_KMC_EXCEPT_CUDA_H
#define KMC_KMC_EXCEPT_CUDA_H

#include <stdexcept>
#include <string>
#include <sstream>

#ifdef __CUDACC__
#include <cuda.h>
#else
typedef unsigned int cudaError_t;
#endif

namespace KmcExcept
{
	class CUDAError : public std::exception
	{
		cudaError_t m_cudaError;
		std::string m_what;

		public:

			CUDAError(cudaError_t cudaError, const char* file, int line, const char* moreWhat = 0);
			
		virtual const char* what() const noexcept {return m_what.c_str();}
#ifdef __CUDACC__
		cudaError_t cudaError() const {return m_cudaError;}
#endif
	};
}
#endif
