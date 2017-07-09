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

#ifndef KPZ_THREAD_LAYOUT_H
#define KPZ_THREAD_LAYOUT_H

#include "systemSize.h"

#include <iostream>

namespace Kpz
{
class ThreadLayout
{
	public:
	static const int WORK_DIM = 2;
	static const int L_MCS_DIV = 1;

	protected:
	size_t m_localWorkSize[WORK_DIM];
	size_t m_globalWorkSize[WORK_DIM];
	int m_blocks;

	int m_lThreadsX, m_lThreadsY;
	int m_lThreadCellDimXW, m_lThreadCellDimYW;

	public:
#ifdef KPZ_INNER_DC
	static const int BASE_LMEM_USEAGE_W = 1;
#else
	static const int BASE_LMEM_USEAGE_W = 4;
#endif

		ThreadLayout()
			: m_lThreadsX(0), m_lThreadsY(0), m_lThreadCellDimXW(0), m_lThreadCellDimYW(0), m_blocks(-1) {}
		ThreadLayout(int ltx, int lty, int ltcdx = 1, int ltcdy = 1)
			: m_lThreadsX(ltx), m_lThreadsY(lty), m_lThreadCellDimXW(ltcdx), m_lThreadCellDimYW(ltcdy) {}

	inline int memPerThreadW() const {return 1<<(m_lThreadCellDimXW+m_lThreadCellDimYW);}
	inline int lmemUsageW(int base = BASE_LMEM_USEAGE_W, int perThreadOverhead = 0) const {
		return base + ((memPerThreadW() + perThreadOverhead)<<(m_lThreadCellDimXW + m_lThreadCellDimYW));
	}
	inline const size_t* localWorkSize() const {return m_localWorkSize;}
	inline const size_t* globalWorkSize() const {return m_globalWorkSize;}
	inline int blocks() const {return m_blocks;}
	inline int lBlockDimXW () const {return m_lThreadsX+m_lThreadCellDimXW;}
	inline int lBlockDimYW () const {return m_lThreadsY+m_lThreadCellDimYW;}
	inline int blockDimXW () const {return 1<<(m_lThreadsX+m_lThreadCellDimXW);}
	inline int blockDimYW () const {return 1<<(m_lThreadsY+m_lThreadCellDimYW);}
	inline int blockData () const {return 1<<(lBlockDimXW() + lBlockDimYW());}
	inline int getRandomNumberCount () const {
		return m_globalWorkSize[0]*m_globalWorkSize[1];
	}

	std::ostream& printConst(std::ostream& o) const;
	void setOptimal(int maxThreads, int maxLmem, int blocks, SystemSize& size);

	template<class LocalLayout>
	static ThreadLayout fromLocalLayout() {
		return ThreadLayout(LocalLayout::L_THREADS_X, LocalLayout::L_THREADS_Y
				, LocalLayout::L_THREAD_CELL_DIM_X_W, LocalLayout::L_THREAD_CELL_DIM_Y);
	}
};
}

#endif
