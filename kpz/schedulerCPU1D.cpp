/***************************************************************************
*   Copyright 2013 Jeffrey Kelling <j.kelling@hzdr.de>
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

#include "schedulerCPU1D.h"

#include "kpzConst.h"
#include "systemSize.h"
#include "lattice.h"

#include <kmcRandom.h>
#include <benchmarking.h>

#include <iomanip>
#include <iostream>
#include <cstring>

void Kpz::SchedulerCPU1D::mcs(int n)
{
	for(; n > 0; --n)
	{
		for(int y = 0; y < m_size.dimY(); ++y)
		for(int a = 0; a < m_size.dimX(); ++a)
		{
			const int x = (int)(dsfmt_genrand_close_open(&m_dsfmt)*m_size.dimX());

			Slope c = m_size.slopeX(x, y);
			Slope r = m_size.slopeX((x+1)&m_sizeCache.mDimX(), y);

			if((!c(m_stash)) && r(m_stash))
			{
				c.flip(m_stash);
				r.flip(m_stash);
			}
		}
		++m_deviceMcs;
	}
}

void Kpz::SchedulerCPU1D::mcsPQ(int n)
{
	for(; n > 0; --n)
	{
		for(int y = 0; y < m_size.dimY(); ++y)
		for(int a = 0; a < m_size.dimX(); ++a)
		{
			const int x = (int)(dsfmt_genrand_close_open(&m_dsfmt)*m_size.dimX());

			Slope c = m_size.slopeX(x, y);
			Slope r = m_size.slopeX((x+1)&m_sizeCache.mDimX(), y);

			const int config = c(m_stash) | (r(m_stash)<<1);

			if(config == 0b10)
			{ // p
				if(dsfmt_genrand_close_open(&m_dsfmt) >= m_disorderP[0])
					continue;
			}
			else if(config == 0b01)
			{ // q
				if(dsfmt_genrand_close_open(&m_dsfmt) >= m_disorderQ[0])
					continue;
			}
			else continue;
			c.flip(m_stash);
			r.flip(m_stash);
		}
		++m_deviceMcs;
	}
}

void Kpz::SchedulerCPU1D::mcsDisorder2(int n)
{
	if(!m_disorder)
		return;
	for(; n > 0; --n)
	{
		for(int y = 0; y < m_size.dimY(); ++y)
		for(int a = 0; a < m_size.dimX(); ++a)
		{
			const int x = (int)(dsfmt_genrand_close_open(&m_dsfmt)*m_size.dimX());

			Slope c = m_size.slopeX(x, y);
			Slope r = m_size.slopeX((x+1)&m_sizeCache.mDimX(), y);

			const int config = c(m_stash) | (r(m_stash)<<1);

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
		}
		++m_deviceMcs;
	}
}

void Kpz::SchedulerCPU1D::mcsSyncRng(int n)
{
	for(; n > 0; --n)
	{
		for(int y = 0; y < m_size.dimY(); ++y)
		for(int a = 0; a < m_size.dimX(); ++a)
		{
			RngPackage rnd(m_dsfmt);
			const int x = (int)(rnd[0]*m_size.dimX());

			Slope c = m_size.slopeX(x, y);
			Slope r = m_size.slopeX((x+1)&m_sizeCache.mDimX(), y);

			if((!c(m_stash)) && r(m_stash))
			{
				c.flip(m_stash);
				r.flip(m_stash);
			}
		}
		++m_deviceMcs;
	}
}

void Kpz::SchedulerCPU1D::mcsPQSyncRng(int n)
{
	for(; n > 0; --n)
	{
		for(int y = 0; y < m_size.dimY(); ++y)
		for(int a = 0; a < m_size.dimX(); ++a)
		{
			RngPackage rnd(m_dsfmt);
			const int x = (int)(rnd[0]*m_size.dimX());

			Slope c = m_size.slopeX(x, y);
			Slope r = m_size.slopeX((x+1)&m_sizeCache.mDimX(), y);

			const int config = c(m_stash) | (r(m_stash)<<1);

			if(config == 0b10)
			{ // p
				if(rnd[2] >= m_disorderP[0])
					continue;
			}
			else if(config == 0b01)
			{ // q
				if(rnd[2] >= m_disorderQ[0])
					continue;
			}
			else continue;
			c.flip(m_stash);
			r.flip(m_stash);
		}
		++m_deviceMcs;
	}
}

void Kpz::SchedulerCPU1D::mcsDisorder2SyncRng(int n)
{
	if(!m_disorder)
		return;
	for(; n > 0; --n)
	{
		for(int y = 0; y < m_size.dimY(); ++y)
		for(int a = 0; a < m_size.dimX(); ++a)
		{
			RngPackage rnd(m_dsfmt);
			const int x = (int)(rnd[0]*m_size.dimX());

			Slope c = m_size.slopeX(x, y);
			Slope r = m_size.slopeX((x+1)&m_sizeCache.mDimX(), y);

			const int config = c(m_stash) | (r(m_stash)<<1);

			if(config == 0b10)
			{ // p
				if(rnd[2] >= m_disorderP[c(m_disorder)])
					continue;
			}
			else if(config == 0b01)
			{ // q
				if(rnd[2] >= m_disorderQ[c(m_disorder)])
					continue;
			}
			else continue;
			c.flip(m_stash);
			r.flip(m_stash);
		}
		++m_deviceMcs;
	}
}

void Kpz::SchedulerCPU1D::correlateTag()
{
	join();
	pthread_create(&m_roughnessThread, NULL, &threadCorrelateLinesTag, (void*) this);
}
