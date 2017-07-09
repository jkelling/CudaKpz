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

#include "LocalLayout.h"

#include <iostream>

namespace Kpz
{
	template<class LocalLayout>
	std::ostream& printLocalLayoutBasics(std::ostream& out)
	{
		out << "\tThreads (x,y): " << LocalLayout::THREADS_X << ' ' << LocalLayout::THREADS_Y
			<< " total: " << LocalLayout::THREADS
			<< "\n\tThreadCellDimW (x,y): " << LocalLayout::THREAD_CELL_DIM_X_W << ' ' << LocalLayout::THREAD_CELL_DIM_Y_W
			<< " size: " << (1<<(LocalLayout::L_THREAD_CELL_DIM_X_W+LocalLayout::L_THREAD_CELL_DIM_Y_W))
			<< "\n\tSitesPerThreadCell: " << LocalLayout::SITES_PER_THREAD_CELL
			<< "\n\tL_MCS_DIV: " << LocalLayout::L_MCS_DIV << " , updatesPerThreadCell: " << LocalLayout::UPDATES_PER_THREAD_CELL
			<< "\n\tBlockDimW: (x,y): " << LocalLayout::BLOCK_DIM_X_W << ' ' << LocalLayout::BLOCK_DIM_Y_W
			<< " sizeW: " << LocalLayout::BLOCK_DATA
			;
		return out;
	}

	template<int lThreadsX, int lThreadsY, int lMCSdiv, int lThreadCellDimXW, int lThreadCellDimYW>
	std::ostream& LocalLayout<lThreadsX, lThreadsY, lMCSdiv, lThreadCellDimXW, lThreadCellDimYW>::
	print(std::ostream& out)
	{
		out << "LocalLayout:\n";
		return printLocalLayoutBasics<LocalLayout<lThreadsX, lThreadsY, lMCSdiv, lThreadCellDimXW, lThreadCellDimYW> >
			(out) << std::endl;
	}
	
	template<int lThreadsX, int lThreadsY, int lMCSdiv, int lThreadCellDimXW, int lThreadCellDimYW>
	std::ostream& LocalLayoutDB4<lThreadsX, lThreadsY, lMCSdiv, lThreadCellDimXW, lThreadCellDimYW>::
	print(std::ostream& out)
	{
		out << "LocalLayoutDB4:\n";
		printLocalLayoutBasics<LocalLayoutDB4<lThreadsX, lThreadsY, lMCSdiv, lThreadCellDimXW, lThreadCellDimYW> >(out);
		out << "\nBlockBorderWidth : " << BLOCK_BORDER_WIDTH << std::endl;
		return out;
	}

	template<int lThreadsX, int lThreadsY, int lMCSdiv, int lThreadCellDimXW, int lThreadCellDimYW>
	std::ostream& LocalLayoutDT<lThreadsX, lThreadsY, lMCSdiv, lThreadCellDimXW, lThreadCellDimYW>::
	print(std::ostream& out)
	{
		out << "LocalLayoutDT:\n";
		return printLocalLayoutBasics<LocalLayoutDT<lThreadsX, lThreadsY, lMCSdiv, lThreadCellDimXW, lThreadCellDimYW> >
			(out) << std::endl;
	}
}
