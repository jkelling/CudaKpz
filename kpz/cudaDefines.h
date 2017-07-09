/***************************************************************************
*   Copyright 2011 - 2014 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef KPZ_CUDA_DEFINES_H
#define KPZ_CUDA_DEFINES_H

#define GPU_BLOCK_ID_X blockIdx.x
#define GPU_THREAD_ID_X threadIdx.x
#define GPU_THREAD_ID_Y threadIdx.y
#define GPU_BLOCK_DIM_X blockDim.x
#define GPU_BLOCK_DIM_Y blockDim.y
#define GPU_GRID_DIM_X gridDim.x
#define GPU_SHARED_MEM __shared__
#define GPU_LOCAL_SYNC __syncthreads();
#define GPU_CAST(space,type,val) ((type)(val))

#endif
