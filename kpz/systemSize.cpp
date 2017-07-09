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

#include "systemSize.h"

#include "threadLayout.h"

#include <iostream>
#include <utility>

void Kpz::SystemSize::adjust(const ThreadLayout& layout)
{
	m_lBlocksX = m_lDimXW - layout.lBlockDimXW();
	const int lBlocksY = m_lDimYW - layout.lBlockDimYW();
	if(m_lBlocksX < 0 || lBlocksY < 0)
		m_mBlocksGPU = -1;
	else
		m_mBlocksGPU = (1<<(m_lBlocksX+lBlocksY))-1;
}

bool Kpz::SystemSizeBase::set(int lx, int ly)
{
	//TODO: check for minumum size
	if(lx > 30 || ly > 30)
	{
		std::cerr << "Systemdimensions are to large.\n";
		return false;
	}
	m_lDimX = lx;
	m_lDimY = ly;
	m_dimX = 1<<lx;
	m_dimY = 1<<ly;
	m_lDimXW = lx-L_BASIC_CELL_DIM_X;
	m_lDimYW = ly-L_BASIC_CELL_DIM_Y;
	m_mDimXW = (1<<m_lDimXW)-1;
	m_mDimYW = (1<<m_lDimYW)-1;

	m_lBlocksX = 0;
	m_mBlocksGPU = 0;
}

void Kpz::SystemSizeCPU::cpuSet(const SystemSize& size)
{
	m_mDimX = size.dimX()-1;
	m_mDimY = size.dimY()-1;
}

unsigned int* Kpz::SystemSize::init() const
{
	return new unsigned int[sizeW()];
}

std::ostream& Kpz::operator<<(std::ostream& o, const Kpz::SystemSize& s)
{
	o << "SystemSize = " << s.dimX() << '*' << s.dimY()
	<< " = " << s.size()
	<< "\nSystem Data = " << s.sizeW() << " quad = " << (s.sizeW()<<2) << " byte"
	<< "\nDimensions: " << s.dimXW() << 'x' << s.dimYW() << "; Row = " << s.dimXW()
	<< "\nNumber of Blocks = " << s.blocksGPU()
	<< "\nWords per row (x) = " << s.dimXW() << ", Blocks: " << (1 << s.lBlocksX())
	<< "\nNumber of rows (y) = " << s.dimYW() << ", Blocks: " << (s.blocksGPU()>>s.lBlocksX())
	<< std::endl;
}

std::ostream& Kpz::SystemSize::printConst(std::ostream& o) const
{
	return o << "#define SIZE_M_BLOCKS_GPU " << m_mBlocksGPU
		<< "\n#define SIZE_L_BLOCKS_X " << m_lBlocksX
		<< "\n#define SIZE_L_DIM_X " << m_lDimX
		<< "\n#define SIZE_L_DIM_Y " << m_lDimY
		<< "\n#define SIZE_DIM_X " << m_dimX
		<< "\n#define SIZE_DIM_Y " << m_dimY
		<< "\n#define SIZE_L_DIM_X_W " << m_lDimXW
		<< "\n#define SIZE_L_DIM_Y_W " << m_lDimYW
		<< "\n#define SIZE_M_DIM_X_W " << m_mDimXW
		<< "\n#define SIZE_M_DIM_Y_W " << m_mDimYW
		<< "\n\n";
}

void Kpz::SystemSize::adjustDT()
{
	m_lBlocksX-=1;
	m_mBlocksGPU>>=2;
	if(m_mBlocksGPU == 0 || m_lBlocksX < 0)
		m_mBlocksGPU = -1;
}

#include <splash/splash.h>

splash::Dimensions Kpz::SystemSize::splashDimensions() const
{
	return std::move(splash::Dimensions(dimXW(), dimYW(), 1));
}
