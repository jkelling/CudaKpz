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

#include "threadLayout.h"

#include <algorithm>

void Kpz::ThreadLayout::setOptimal(int maxThreads, int maxLmem, int blocks, SystemSize& size)
{
	switch(maxThreads) //TODO
	{
		case 1024:
			m_lThreadsX = m_lThreadsY = 5;
			m_lThreadCellDimXW = 2;
			m_lThreadCellDimYW = 1;
			break;
		case 512:
			m_lThreadsX = 5;
			m_lThreadsY = 4;
			m_lThreadCellDimXW = m_lThreadCellDimYW = 1;
			break;
		case 256:
			m_lThreadsX = m_lThreadsY = 4;
			m_lThreadCellDimXW = m_lThreadCellDimYW = 2;
			break;
	}
	size.adjust(*this);
	m_blocks = std::min(blocks, size.blocksGPU());

	m_localWorkSize[0] = 1 << m_lThreadsX;
	m_localWorkSize[1] = 1 << m_lThreadsY;
	m_globalWorkSize[0] = m_localWorkSize[0] * m_blocks;
	m_globalWorkSize[1] = m_localWorkSize[1];
}

std::ostream& Kpz::ThreadLayout::printConst(std::ostream& o) const
{
	const int THREADS_X = 1 << m_lThreadsX;
	const int THREADS_Y = 1 << m_lThreadsY;
	const int L_THREADS = m_lThreadsX + m_lThreadsY;
	const int THREADS = 1<<L_THREADS;

	const int L_BASIC_CELL_DIM_X = 2; // double word, 2 bits per site
	const int BASIC_CELL_DIM_X = 1<<L_BASIC_CELL_DIM_X;
	const int L_BASIC_CELL_DIM_Y = 2;
	const int BASIC_CELL_DIM_Y = 1<<L_BASIC_CELL_DIM_Y;
	const int L_THREAD_CELL_DIM_X = (L_BASIC_CELL_DIM_X + m_lThreadCellDimXW);
	const int L_THREAD_CELL_DIM_Y = (L_BASIC_CELL_DIM_Y + m_lThreadCellDimYW);
	const int THREAD_CELL_DIM_X = 1<<(L_THREAD_CELL_DIM_X);
	const int THREAD_CELL_DIM_Y = 1<<(L_THREAD_CELL_DIM_Y);
	const int THREAD_CELL_DIM_X_W = 1<<(m_lThreadCellDimXW);
	const int THREAD_CELL_DIM_Y_W = 1<<(m_lThreadCellDimYW);
	const int DEAD_BORDER_DIM = 1;
	const int ACT_THREAD_CELL_DIM_X = THREAD_CELL_DIM_X - DEAD_BORDER_DIM;
	const int ACT_THREAD_CELL_DIM_Y = THREAD_CELL_DIM_Y - DEAD_BORDER_DIM;
	const int L_BLOCK_DIM_X_W = lBlockDimXW();
	const int L_BLOCK_DIM_X = L_BLOCK_DIM_X_W + L_BASIC_CELL_DIM_X;
	const int BLOCK_DIM_X_W = 1<<L_BLOCK_DIM_X_W;
	const int M_BLOCK_DIM_X = (1<<L_BLOCK_DIM_X)-1;
	const int M_BLOCK_DIM_X_W = BLOCK_DIM_X_W-1;
	const int L_BLOCK_DIM_Y_W = lBlockDimYW();
	const int L_BLOCK_DIM_Y = L_BLOCK_DIM_Y_W + L_BASIC_CELL_DIM_Y;
	const int BLOCK_DIM_Y_W = 1<<L_BLOCK_DIM_Y_W;
	const int BLOCK_BORDER_Y = (1<<L_BLOCK_DIM_Y)-1;
	const int BLOCK_BORDER_X = (1<<L_BLOCK_DIM_X)-1;
	const int M_BLOCK_DIM_Y = (1<<L_BLOCK_DIM_Y)-1;
	const int M_BLOCK_DIM_Y_W = BLOCK_DIM_Y_W-1;
	const int L_BLOCK_DATA = L_BLOCK_DIM_X_W+L_BLOCK_DIM_Y_W;
	const int BLOCK_DATA = 1<<L_BLOCK_DATA;
	const int BLOCK_SIZE_W = BLOCK_DATA;
	const int M_BLOCK_DATA = BLOCK_DATA-1;
	
	const int SITES_PER_THREAD_CELL = THREAD_CELL_DIM_X*THREAD_CELL_DIM_Y;
	const int ACTIVE_SITES_PER_THREAD_CELL = (THREAD_CELL_DIM_X-1)*(THREAD_CELL_DIM_Y-1);
	const int UPDATES_PER_THREAD_CELL = SITES_PER_THREAD_CELL>>L_MCS_DIV;

	//KPZ_INNER_DC
	const int L_SUB_THREAD_CELL_DIM_X = (L_THREAD_CELL_DIM_X-1);
	const int L_SUB_THREAD_CELL_DIM_Y = (L_THREAD_CELL_DIM_Y-1);
	const int M_SUB_THREAD_CELL_DIM_X = (1 << L_SUB_THREAD_CELL_DIM_X)-1;
	const int M_SUB_THREAD_CELL_DIM_Y = (1 << L_SUB_THREAD_CELL_DIM_Y)-1;

	return o << "#define L_THREADS_X " << m_lThreadsX
		<< "\n#define L_THREADS_Y " << m_lThreadsY
		<< "\n#define THREADS_X " << THREADS_X
		<< "\n#define THREADS_Y " << THREADS_Y
		<< "\n#define L_THREADS " << L_THREADS
		<< "\n#define L_BASIC_CELL_DIM_X " << L_BASIC_CELL_DIM_X
		<< "\n#define L_BASIC_CELL_DIM_Y " << L_BASIC_CELL_DIM_Y
		<< "\n#define BASIC_CELL_DIM_X " << BASIC_CELL_DIM_X
		<< "\n#define BASIC_CELL_DIM_Y " << BASIC_CELL_DIM_Y
		<< "\n#define L_THREAD_CELL_DIM_X " << L_THREAD_CELL_DIM_X
		<< "\n#define L_THREAD_CELL_DIM_Y " << L_THREAD_CELL_DIM_Y
		<< "\n#define L_THREAD_CELL_DIM_X_W " << m_lThreadCellDimXW
		<< "\n#define L_THREAD_CELL_DIM_Y_W " << m_lThreadCellDimYW
		<< "\n#define THREAD_CELL_DIM_X " << THREAD_CELL_DIM_X
		<< "\n#define THREAD_CELL_DIM_Y " << THREAD_CELL_DIM_Y
		<< "\n#define THREAD_CELL_DIM_X_W " << THREAD_CELL_DIM_X_W
		<< "\n#define THREAD_CELL_DIM_Y_W " << THREAD_CELL_DIM_Y_W
		<< "\n#define DEAD_BORDER_DIM " << DEAD_BORDER_DIM
		<< "\n#define ACT_THREAD_CELL_DIM_X " << ACT_THREAD_CELL_DIM_X
		<< "\n#define ACT_THREAD_CELL_DIM_Y " << ACT_THREAD_CELL_DIM_Y
		<< "\n#define L_BLOCK_DIM_X " << L_BLOCK_DIM_X
		<< "\n#define L_BLOCK_DIM_Y " << L_BLOCK_DIM_Y
		<< "\n#define L_BLOCK_DIM_X_W " << L_BLOCK_DIM_X_W
		<< "\n#define L_BLOCK_DIM_Y_W " << L_BLOCK_DIM_Y_W
		<< "\n#define BLOCK_DIM_X_W " << BLOCK_DIM_X_W
		<< "\n#define BLOCK_DIM_Y_W " << BLOCK_DIM_Y_W
		<< "\n#define M_BLOCK_DIM_X " << M_BLOCK_DIM_X
		<< "\n#define M_BLOCK_DIM_X_W " << M_BLOCK_DIM_X_W
		<< "\n#define M_BLOCK_DIM_Y " << M_BLOCK_DIM_Y
		<< "\n#define M_BLOCK_DIM_Y_W " << M_BLOCK_DIM_Y_W
		<< "\n#define BLOCK_BORDER_X " << BLOCK_BORDER_X
		<< "\n#define BLOCK_BORDER_Y " << BLOCK_BORDER_Y
		<< "\n#define L_BLOCK_DATA " << L_BLOCK_DATA
		<< "\n#define BLOCK_DATA " << BLOCK_DATA
		<< "\n#define M_BLOCK_DATA " << M_BLOCK_DATA
		<< "\n#define SITES_PER_THREAD_CELL " << SITES_PER_THREAD_CELL
		<< "\n#define ACTIVE_SITES_PER_THREAD_CELL " << ACTIVE_SITES_PER_THREAD_CELL
		<< "\n#define L_MCS_DIV " << L_MCS_DIV
		<< "\n#define UPDATES_PER_THREAD_CELL " << UPDATES_PER_THREAD_CELL
		<< "\n#define L_SUB_THREAD_CELL_DIM_X " << L_SUB_THREAD_CELL_DIM_X
		<< "\n#define L_SUB_THREAD_CELL_DIM_Y " << L_SUB_THREAD_CELL_DIM_Y
		<< "\n#define M_SUB_THREAD_CELL_DIM_X " << M_SUB_THREAD_CELL_DIM_X
		<< "\n#define M_SUB_THREAD_CELL_DIM_Y " << M_SUB_THREAD_CELL_DIM_Y
		<< "\n\n";
}
