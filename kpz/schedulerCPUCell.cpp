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

#include "schedulerCPUCell.h"

#include "kpzConst.h"
#include "systemSize.h"
#include "LocalLayoutIO.h"

#include <kmcRandom.h>

#include <iomanip>
#include <iostream>
#include <list>

template<class LocalLayout>
void Kpz::SchedulerCPUCell<LocalLayout>::mcs(int n)
{
	for(; n > 0; --n)
	{
		const int xw = (int)(dsfmt_genrand_close_open(&m_dsfmt)*MyLocalLayout::BLOCK_DIM_X_W);
		const int yw = (int)(dsfmt_genrand_close_open(&m_dsfmt)*MyLocalLayout::BLOCK_DIM_Y_W);
		for(int blocks = 0; blocks < m_size.sizeW()>>MyLocalLayout::L_BLOCK_DATA; ++blocks)
		{
			updateThreadCell(blocks, xw, yw);
		}
		++m_deviceMcs;
	}
}

template<class LocalLayout>
void Kpz::SchedulerCPUCell<LocalLayout>::updateThreadCell(int blocks, int xw, int yw)
{
	int lSystem[MyLocalLayout::BLOCK_DIM_X_W*MyLocalLayout::BLOCK_DIM_Y_W];

	//load
	xw = ((xw + (blocks << MyLocalLayout::L_BLOCK_DIM_X_W)) & m_size.mDimXW());
	yw = ((yw + ((blocks>>m_size.lBlocksX()) << MyLocalLayout::L_BLOCK_DIM_Y_W)) & m_size.mDimYW());
	for (int a = 0, y = yw; a<MyLocalLayout::BLOCK_DIM_Y_W; ++a, (++y)&=m_size.mDimYW())
	{
		for (int s = 0, x = (xw)&m_size.mDimXW()
			; s < MyLocalLayout::BLOCK_DIM_X_W; s+=1, (x+=1)&=m_size.mDimXW())
		{
			lSystem[(a<<MyLocalLayout::L_BLOCK_DIM_X_W) + s] = m_stash[(y<<m_size.lDimXW()) + x];
		}
	}

	updateThreadCell(lSystem);

	//save
	for (int a = 0, y = yw; a<MyLocalLayout::BLOCK_DIM_Y_W; ++a, (++y)&=m_size.mDimYW())
	{
		for (int s = 0, x = (xw)&m_size.mDimXW()
			; s < MyLocalLayout::BLOCK_DIM_X_W; s+=1, (x+=1)&=m_size.mDimXW())
		{
			m_stash[(y<<m_size.lDimXW()) + x] = lSystem[(a<<MyLocalLayout::L_BLOCK_DIM_X_W) + s];
		}
	}
}

namespace Kpz
{
	struct Delay {
		int xIndex, xShift, yIndex, yShift, thisIndex, thisShift;

			Delay (int xIndex, int xShift, int yIndex, int yShift, int thisIndex, int thisShift)
				: xIndex(xIndex), xShift(xShift), yIndex(yIndex), yShift(yShift), thisIndex(thisIndex), thisShift(thisShift)
			{}
	};
}

#define WAIT_IF_DEAD

template<class LocalLayout>
void Kpz::SchedulerCPUCell<LocalLayout>::updateThreadCell(int* system)
{

	#ifdef WAIT_IF_DEAD
		std::list<Delay> wait;
		std::list<Delay> waitXY;
	#endif
	int ref[2];
	for(int a = 0; a < MyLocalLayout::UPDATES_PER_THREAD_CELL; ++a) //NOTE number of acive sites is lower
	{
		//using lowest 2 * 12 bit for offest
		for(int tx = 0; tx < MyLocalLayout::THREADS_X; ++tx)
			for(int ty = 0; ty < MyLocalLayout::THREADS_Y; ++ty)
			{
				const int dispX = tx << MyLocalLayout::L_THREAD_CELL_DIM_X_W << L_BASIC_CELL_DIM_X;
				const int dispY = ty << MyLocalLayout::L_THREAD_CELL_DIM_Y_W << L_BASIC_CELL_DIM_Y;
				int x = dsfmt_genrand_close_open(&m_dsfmt)*KMC_LCG_RAND_SUP;
				int y = (x>>MyLocalLayout::O_WIDTH)&MyLocalLayout::O_MASK;
				x &= MyLocalLayout::O_MASK;
				if(!(tx | ty)) // thread 0 picks reference
				{
					ref[0] = ((x>>L_BASIC_CELL_DIM_X)&MyLocalLayout::M_BLOCK_DIM_X_W)<<L_BASIC_CELL_DIM_X;
					ref[1] = ((y>>L_BASIC_CELL_DIM_Y)&MyLocalLayout::M_BLOCK_DIM_Y_W)<<L_BASIC_CELL_DIM_Y;
				}
				#ifdef WAIT_IF_DEAD
				const int lx = (x&MyLocalLayout::M_THREAD_CELL_DIM_X);
				const int ly = (y&MyLocalLayout::M_THREAD_CELL_DIM_Y);
				x = (lx + ref[0] + dispX) & MyLocalLayout::M_BLOCK_DIM_X;
				y = (ly + ref[1] + dispY) & MyLocalLayout::M_BLOCK_DIM_Y;
				#else
				x &= MyLocalLayout::O_MASK;
				x = ((int)((float)x/MyLocalLayout::O_RND_SUP*ACT_THREAD_CELL_DIM_X) + ref[0] + dispX) & MyLocalLayout::M_BLOCK_DIM_X;
				y = ((int)((float)y/MyLocalLayout::O_RND_SUP*ACT_THREAD_CELL_DIM_Y) + ref[1] + dispY) & MyLocalLayout::M_BLOCK_DIM_Y;
				#endif
				//std::cerr << "\n#################################\nAttempting to update ( " << x << " , " << y << " )\n";

				if((x == MyLocalLayout::BLOCK_BORDER_X) || (y == MyLocalLayout::BLOCK_BORDER_Y))
					continue;

				const int thisIndex = MyLocalLayout::index(x, y);
				const int thisShift = MyLocalLayout::shift(x, y);

				//std::cerr << "condition 1: " << ((system[thisIndex]>>thisShift)&3) << " == " << 3;
				// p = 1, q = 0 for now
				if(((system[thisIndex]>>thisShift)&3) != 0) 
					continue; // no p rocess possible

				const int xIndex = MyLocalLayout::index((x+1), y);
				const int xShift = MyLocalLayout::shift((x+1), y);
				const int yIndex = MyLocalLayout::index(x, (y+1));
				const int yShift = MyLocalLayout::shift(x, (y+1)) | 1; // select slope in y direction

				//std::cerr << " -- satisfied\ncondition 2: !(" << ((system[xIndex]>>xShift)&1) << 
					//" | " << ((system[yIndex]>>yShift)&1) << ')';

				if(!(((system[xIndex]>>xShift)&(system[yIndex]>>yShift))&1))
				//if((((system[xIndex]>>xShift)|(system[yIndex]>>yShift))&1))
					continue;

				#ifdef WAIT_IF_DEAD
				if((lx == MyLocalLayout::M_THREAD_CELL_DIM_X) && (ly == MyLocalLayout::M_THREAD_CELL_DIM_Y))
				{
					waitXY.push_back(Delay(xIndex, xShift, yIndex, yShift, thisIndex, thisShift));
					continue;
				}
				else if ((lx == MyLocalLayout::M_THREAD_CELL_DIM_X) || (ly == MyLocalLayout::M_THREAD_CELL_DIM_Y))
				{
					wait.push_back(Delay(xIndex, xShift, yIndex, yShift, thisIndex, thisShift));
					continue;
				}
				#endif
				//std::cerr << " -- satisfied\n";

				system[thisIndex] ^= 3<<thisShift;
				system[xIndex] ^= 1<<xShift;
				system[yIndex] ^= 1<<yShift;
			}
		#ifdef WAIT_IF_DEAD
		for(std::list<Delay>::iterator a = wait.begin(); wait.size(); a=wait.erase(a))
		{
			system[a->thisIndex] ^= 3<<a->thisShift;
			system[a->xIndex] ^= 1<<a->xShift;
			system[a->yIndex] ^= 1<<a->yShift;
		}
		for(std::list<Delay>::iterator a = waitXY.begin(); waitXY.size(); a=waitXY.erase(a))
		{
			system[a->thisIndex] ^= 3<<a->thisShift;
			system[a->xIndex] ^= 1<<a->xShift;
			system[a->yIndex] ^= 1<<a->yShift;
		}
		#endif
	}
	//std::cerr << '\n';
}

template<class LocalLayout>
void Kpz::SchedulerCPUCell<LocalLayout>::init()
{
	std::cout << "SchedulerCPUCell (CPU, DB ur, inner DB)\n";
	MyLocalLayout::print(std::cout);
	std::cerr << "[WW][SchedulerCPUCell] This scheduler was not tested"
		"after template<class LocalLayout> refactroization.\n";
}

#include "schedulerCPUCellInstances.h"
