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

#pragma once

#include "kpzConst.h"

#include <iostream>

#ifdef __CUDACC__
#include <cuda.h>
#else
#define __device__
#endif

namespace Kpz
{
	template<int lThreadsX, int lThreadsY, int lMCSdiv, int lThreadCellDimXW, int lThreadCellDimYW>
	class LocalLayout
	{
		public:

		static const int L_THREADS_X = lThreadsX;
		static const int L_THREADS_Y = lThreadsY;
		static const int L_MCS_DIV = lMCSdiv;
		static const int L_THREAD_CELL_DIM_X_W = lThreadCellDimXW;
		static const int L_THREAD_CELL_DIM_Y_W = lThreadCellDimYW;

		static const int L_THREADS = L_THREADS_X + L_THREADS_Y;
		static const int THREADS_X = 1 << L_THREADS_X;
		static const int THREADS_Y = 1 << L_THREADS_Y;
		static const int THREADS = 1<<L_THREADS;
		
		static const int L_BASIC_CELL_DIM_X = Kpz::L_BASIC_CELL_DIM_X; // double word, 2 bits per site
		static const int BASIC_CELL_DIM_X = 1<<L_BASIC_CELL_DIM_X;
		static const int L_BASIC_CELL_DIM_Y = Kpz::L_BASIC_CELL_DIM_Y;
		static const int BASIC_CELL_DIM_Y = 1<<L_BASIC_CELL_DIM_Y;
#ifdef __CUDACC__
		static dim3 threadConfig() {return dim3(THREADS_X, THREADS_Y, 1);}
#endif

		//all dimensions are those of the slope-grid
		static const int L_THREAD_CELL_DIM_X = (Kpz::L_BASIC_CELL_DIM_X + L_THREAD_CELL_DIM_X_W);
		static const int L_THREAD_CELL_DIM_Y = (Kpz::L_BASIC_CELL_DIM_Y + L_THREAD_CELL_DIM_Y_W);
		static const int THREAD_CELL_DIM_X = 1<<(L_THREAD_CELL_DIM_X);
		static const int THREAD_CELL_DIM_Y = 1<<(L_THREAD_CELL_DIM_Y);
		static const int THREAD_CELL_DIM_X_W = 1<<(L_THREAD_CELL_DIM_X_W);
		static const int THREAD_CELL_DIM_Y_W = 1<<(L_THREAD_CELL_DIM_Y_W);
		static const int SITES_PER_THREAD_CELL = THREAD_CELL_DIM_X*THREAD_CELL_DIM_Y;

		static const int L_BLOCK_DIM_X_W = L_THREADS_X+L_THREAD_CELL_DIM_X_W;
		static const int L_BLOCK_DIM_X = L_BLOCK_DIM_X_W + L_BASIC_CELL_DIM_X;
		static const int BLOCK_DIM_X_W = 1<<L_BLOCK_DIM_X_W;
		static const int M_BLOCK_DIM_X = (1<<L_BLOCK_DIM_X)-1;
		static const int M_BLOCK_DIM_X_W = BLOCK_DIM_X_W-1;
		static const int L_BLOCK_DIM_Y_W = L_THREADS_Y+L_THREAD_CELL_DIM_Y_W;
		static const int L_BLOCK_DIM_Y = L_BLOCK_DIM_Y_W + L_BASIC_CELL_DIM_Y;
		static const int BLOCK_DIM_Y_W = 1<<L_BLOCK_DIM_Y_W;

		static const int M_BLOCK_DIM_Y = (1<<L_BLOCK_DIM_Y)-1;
		static const int M_BLOCK_DIM_Y_W = BLOCK_DIM_Y_W-1;
		static const int L_BLOCK_DATA = L_BLOCK_DIM_X_W+L_BLOCK_DIM_Y_W;
		/*! Actual size of block in smem, including borders. May be overridden by derived classes. */
		static const int BLOCK_DATA = 1<<L_BLOCK_DATA;
		/*! Size of active block exculding borders. No need to override in derived classes. */
		static const int BLOCK_SIZE_W = BLOCK_DATA;
		static const int M_BLOCK_SIZE_W = BLOCK_SIZE_W-1;

		static const int UPDATES_PER_THREAD_CELL = SITES_PER_THREAD_CELL>>L_MCS_DIV;

		//KPZ_INNER_DC
		static const int L_SUB_THREAD_CELL_DIM_X = (L_THREAD_CELL_DIM_X-1);
		static_assert(L_SUB_THREAD_CELL_DIM_X <= 8, "Sub-thread cells are too large in x direction (max is 256).");
		static const int L_SUB_THREAD_CELL_DIM_Y = (L_THREAD_CELL_DIM_Y-1);
		static_assert(L_SUB_THREAD_CELL_DIM_Y <= 8, "Sub-thread cells are too large in x direction (max is 256).");
		static const int M_SUB_THREAD_CELL_DIM_X = (1 << L_SUB_THREAD_CELL_DIM_X)-1;
		static const int M_SUB_THREAD_CELL_DIM_Y = (1 << L_SUB_THREAD_CELL_DIM_Y)-1;

		// bitselectionmask for coordinatepicking
		static const int O_MASK = 0xFF;
		static constexpr float O_RND_SUP = O_MASK+1;
		static const int O_WIDTH = 8;
		static const int M_THREAD_CELL_DIM_X_W = THREAD_CELL_DIM_X_W-1;
		static const int M_THREAD_CELL_DIM_Y_W = THREAD_CELL_DIM_Y_W-1;
		static const int M_THREAD_CELL_DIM_X = THREAD_CELL_DIM_X-1;
		static const int M_THREAD_CELL_DIM_Y = THREAD_CELL_DIM_Y-1;

		__device__ static  int index(int x, int y)
		{
			return (y>>L_BASIC_CELL_DIM_Y<<L_BLOCK_DIM_X_W) + (x>>L_BASIC_CELL_DIM_X);
		}

		__device__ static int shift(int x, int y)
		{
			return ((x&3)<<1) // select row by x
				| ((y&3)<<3); // select column by y
		}

		static std::ostream& print(std::ostream&);
	};

	template<int lThreadsX, int lThreadsY, int lMCSdiv, int lThreadCellDimXW, int lThreadCellDimYW>
	class LocalLayoutDB4 : public LocalLayout<lThreadsX, lThreadsY, lMCSdiv, lThreadCellDimXW, lThreadCellDimYW>
	{
		typedef LocalLayout<lThreadsX, lThreadsY, lMCSdiv, lThreadCellDimXW, lThreadCellDimYW> Base;
		public:

		static const int BLOCK_BORDER_WIDTH = 5;
		static const int BLOCK_BORDER_Y = (1<<Base::L_BLOCK_DIM_Y)-BLOCK_BORDER_WIDTH;
		static const int BLOCK_BORDER_X = (1<<Base::L_BLOCK_DIM_X)-BLOCK_BORDER_WIDTH;

		static std::ostream& print(std::ostream&);
	};

	template<int lThreadsX, int lThreadsY, int lMCSdiv, int lThreadCellDimXW, int lThreadCellDimYW>
	class LocalLayoutDT : public LocalLayout<lThreadsX, lThreadsY, lMCSdiv, lThreadCellDimXW, lThreadCellDimYW>
	{
		typedef LocalLayout<lThreadsX, lThreadsY, lMCSdiv, lThreadCellDimXW, lThreadCellDimYW> Base;
		public:

		static const int BLOCK_DATA = (Base::BLOCK_DIM_X_W+1)*(Base::BLOCK_DIM_Y_W+1);

		__device__ static int index(int x, int y)
		{
			return (y>>Base::L_BASIC_CELL_DIM_Y)*(Base::BLOCK_DIM_X_W+1) + (x>>Base::L_BASIC_CELL_DIM_X);
		}

		static std::ostream& print(std::ostream&);
	};

	template<class LocalLayoutBase, template<int, int, int, int, int> class LocalLayoutDerived>
	class LocalLayoutProxy : public LocalLayoutDerived<
			LocalLayoutBase::L_THREADS_X, LocalLayoutBase::L_THREADS_Y
			, LocalLayoutBase::L_MCS_DIV
			, LocalLayoutBase::L_THREAD_CELL_DIM_X_W, LocalLayoutBase::L_THREAD_CELL_DIM_Y_W
			>
	{};
}
