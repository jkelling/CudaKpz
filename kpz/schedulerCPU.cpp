/***************************************************************************
*   Copyright 2011 - 2013 Jeffrey Kelling <j.kelling@hzdr.de>
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

#include "schedulerCPU.h"

#include "kpzConst.h"
#include "systemSize.h"
#include "lattice.h"

#include <kmcRandom.h>
#include <benchmarking.h>

#include <iomanip>
#include <iostream>
#include <cstring>

void Kpz::SchedulerCPU::mcs(int n)
{
	mcs(n, [this](int x, int y){pqOneUpdate(x,y);});
}

void Kpz::SchedulerCPU::mcsPQ(int n)
{
	if(m_disorderP[0] == 1. && m_disorderQ[0] == 1.)
		mcs(n, [this](int x, int y){pqOneOneUpdate(x,y);});
	else
		mcs(n, [this](int x, int y){pqUpdate(x,y);});
}

void Kpz::SchedulerCPU::mcsDisorder2(int n)
{
	if(!m_disorder)
		return;
	for(; n > 0; --n)
	{
		for(long long a = 0; a < m_size.size(); ++a)
		{
			const int x = (int)(dsfmt_genrand_close_open(&m_dsfmt)*m_size.dimX());
			const int y = (int)(dsfmt_genrand_close_open(&m_dsfmt)*m_size.dimY());

			LatticePoint c = m_size.latticePoint(x, y);
			Slope r = m_size.slopeX((x+1)&m_sizeCache.mDimX(), y);
			Slope l = m_size.slopeY(x, (y+1)&m_sizeCache.mDimY());

			const int config = c(m_stash) | (r(m_stash)<<2) | (l(m_stash)<<3);

			if(config == 0b1100)
			{ // p
				if(dsfmt_genrand_close_open(&m_dsfmt) >= m_disorderP[c(m_disorder)])
					continue;
			}
			else if(config == 0b0011)
			{ // q
				if(dsfmt_genrand_close_open(&m_dsfmt) >= m_disorderQ[c(m_disorder)])
					continue;
			}
			else continue;
			c.flip(m_stash);
			r.flip(m_stash);
			l.flip(m_stash);
		}
		++m_deviceMcs;
	}
}

void Kpz::SchedulerCPU::mcsSyncRng(int n)
{
	for(; n > 0; --n)
	{
		for(long long a = 0; a < m_size.size(); ++a)
		{
			RngPackage rnd(m_dsfmt);
			const int x = (int)(rnd[0]*m_size.dimX());
			const int y = (int)(rnd[1]*m_size.dimY());

			LatticePoint c = m_size.latticePoint(x, y);
			Slope r = m_size.slopeX((x+1)&m_sizeCache.mDimX(), y);
			Slope l = m_size.slopeY(x, (y+1)&m_sizeCache.mDimY());

			if((!c(m_stash)) && r(m_stash) && l(m_stash))
			{
				c.flip(m_stash);
				r.flip(m_stash);
				l.flip(m_stash);
			}
		}
		++m_deviceMcs;
	}
}

void Kpz::SchedulerCPU::mcsPQSyncRng(int n)
{
	for(; n > 0; --n)
	{
		for(long long a = 0; a < m_size.size(); ++a)
		{
			RngPackage rnd(m_dsfmt);
			const int x = (int)(rnd[0]*m_size.dimX());
			const int y = (int)(rnd[1]*m_size.dimY());

			LatticePoint c = m_size.latticePoint(x, y);
			Slope r = m_size.slopeX((x+1)&m_sizeCache.mDimX(), y);
			Slope l = m_size.slopeY(x, (y+1)&m_sizeCache.mDimY());

			const int config = c(m_stash) | (r(m_stash)<<2) | (l(m_stash)<<3);

			if(config == 0b1100)
			{ // p
				if(rnd[2] >= m_disorderP[0])
					continue;
			}
			else if(config == 0b0011)
			{ // q
				if(rnd[2] >= m_disorderQ[0])
					continue;
			}
			else continue;
			c.flip(m_stash);
			r.flip(m_stash);
			l.flip(m_stash);
		}
		++m_deviceMcs;
	}
}

void Kpz::SchedulerCPU::mcsDisorder2SyncRng(int n)
{
	if(!m_disorder)
		return;
	for(; n > 0; --n)
	{
		for(long long a = 0; a < m_size.size(); ++a)
		{
			RngPackage rnd(m_dsfmt);
			const int x = (int)(rnd[0]*m_size.dimX());
			const int y = (int)(rnd[1]*m_size.dimY());

			LatticePoint c = m_size.latticePoint(x, y);
			Slope r = m_size.slopeX((x+1)&m_sizeCache.mDimX(), y);
			Slope l = m_size.slopeY(x, (y+1)&m_sizeCache.mDimY());

			const int config = c(m_stash) | (r(m_stash)<<2) | (l(m_stash)<<3);

			if(config == 0b1100)
			{ // p
				if(rnd[2] >= m_disorderP[c(m_disorder)])
					continue;
			}
			else if(config == 0b0011)
			{ // q
				if(rnd[2] >= m_disorderQ[c(m_disorder)])
					continue;
			}
			else continue;
			c.flip(m_stash);
			r.flip(m_stash);
			l.flip(m_stash);
		}
		++m_deviceMcs;
	}
}
