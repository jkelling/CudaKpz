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

#include "correlatorLines.h"

#include <cstring>

Kpz::CorrelatorLines::CorrelatorLines(SchedulerService* scheduler)
	: Correlator(scheduler), m_hL(0)
{
}

// void Kpz::CorrelatorLines::set(int mcs, Kpz::Roughness* r)
// {
// 	m_mcs = mcs;
// 	const SystemSize& size = m_scheduler->size();
// 	delete m_snapshot;
// 	delete m_h;
// 	delete m_ch;
// 	delete m_cs;
// 	m_snapshot = new unsigned int[size.sizeW()];
// 	m_h = new double[size.dimY()];
// 	m_ch = new double[size.dimY()];
// 	m_cs = new double[size.dimY()];
// 	memcpy(m_snapshot, m_scheduler->system(), size.sizeW()<<2);
// 	memset(m_h, 0, size.dimY()*sizeof(double));
// 	memset(m_ch, 0, size.dimY()*sizeof(double));
// 	memset(m_cs, 0, size.dimY()*sizeof(double));
// 	for (int yy=0;yy<size.dimY();++yy)
// 	{
// 		double h = 0, rh = 0, rh2 = 0;
// 		for (int xx=1;xx<size.dimX();++xx) 
// 		{
// 			h += (size.slopeX(xx,yy)(m_snapshot)) ? +1 : -1;
// 			rh += h;
// 			rh2 += h*h;
// 		}
// 		if(size.slopeX(0,yy)(m_snapshot))
// 			++h;
// 		else --h;
// 		rh += h;
// 		rh2 += h*h;
// 		rh /= size.dimX();
// 		rh2 /= size.dimX();
// 		m_h[yy] = rh;
// 		if(r)
// 		{
// 			r[yy].h = rh;
// 			r[yy].h2 = rh2;
// 			r[yy].w2 = rh2 - rh*rh;
// 		}
// 	}
// }
Kpz::Roughness Kpz::CorrelatorLines::set(int mcs)
{
	m_mcs = mcs;
	const SystemSize& size = m_scheduler->size();
	delete m_snapshot;
	m_snapshot = new unsigned int[size.sizeW()];
	m_hL = new double[size.dimY()];
	memcpy(m_snapshot, m_scheduler->system(), size.sizeW()<<2);
	Roughness r;
	for (int yy=0;yy<size.dimY();++yy)
	{
		Roughness rl;
		double h = 0;
		for (int xx=1;xx<size.dimX();++xx) 
		{
			h += (size.slopeX(xx,yy)(m_snapshot)) ? +1 : -1;
			rl.addH(h);
		}
		if(size.slopeX(0,yy)(m_snapshot))
			++h;
		else --h;
		rl.addH(h);
		rl.normAndSetW2(size.dimX());
		r += rl;
		m_hL[yy] = rl.h[0];
	}
	r /= size.dimY();
	m_h = r.h[0];
	return r;
}

// void Kpz::CorrelatorLines::correlate(Kpz::Roughness* r)
// {
// 	const SystemSize& size = m_scheduler->size();
// 	const unsigned int* system = m_scheduler->system();
// 	memset(m_ch, 0, size.dimY()*sizeof(double));
// 	memset(m_cs, 0, size.dimY()*sizeof(double));
// 
// 	for (int yy=0;yy<size.dimY();++yy)
// 	{
// 		double hSnap = 0, hSys = 0;
// 		double rh = 0, rh2 = 0;
// 		long long tcs = 0;
// 		for (int xx=1;xx<size.dimX();++xx) 
// 		{
// 			const Slope col = size.slopeX(xx,yy);
// 			const int colSnap = col(m_snapshot);
// 			const int colSys = col(system);
// 			hSys += (colSys) ? +1 : -1;
// 			tcs += ((colSnap) == (colSys)) ? +1 : -1;
// 			hSnap += (colSnap) ? +1 : -1;
// 			rh += hSys;
// 			rh2 += hSys*hSys;
// 			m_ch[yy] += hSnap*hSys;
// 		}
// 		const Slope row = size.slopeX(0,yy);
// 		const int rowSnap = row(m_snapshot);
// 		const int rowSys = row(system);
// 		hSys += (rowSys) ? +1 : -1;
// 		tcs += ((rowSnap) == (rowSys)) ? +1 : -1;
// 		hSnap += (rowSnap) ? +1 : -1;
// 		m_ch[yy] += hSnap*hSys;
// 		rh += hSys;
// 		rh2 += hSys*hSys;
// 		rh /= size.dimX();
// 		rh2 /= size.dimX();
// 		m_ch[yy] = m_ch[yy]/size.dimX() - rh*m_h[yy];
// 		m_cs[yy] = double(tcs)/size.dimX();
// 		if(r)
// 		{
// 			r[yy].h = rh;
// 			r[yy].h2 = rh2;
// 			r[yy].w2 = rh2 - rh*rh;
// 		}
// 	}
// }
Kpz::Roughness Kpz::CorrelatorLines::correlate()
{
	const SystemSize& size = m_scheduler->size();
	const unsigned int* system = m_scheduler->system();
	Roughness r;
	m_ch = m_cs = 0;

	for (int yy=0;yy<size.dimY();++yy)
	{
		double hSnap = 0, hSys = 0;
		Roughness rl;
		double tch = 0;
		long long tcs = 0;
		for (int xx=1;xx<size.dimX();++xx) 
		{
			const Slope col = size.slopeX(xx,yy);
			const int colSnap = col(m_snapshot);
			const int colSys = col(system);
			hSys += (colSys) ? +1 : -1;
			tcs += ((colSnap) == (colSys)) ? +1 : -1;
			hSnap += (colSnap) ? +1 : -1;
			rl.addH(hSys);
			tch += hSnap*hSys;
		}
		const Slope row = size.slopeX(0,yy);
		const int rowSnap = row(m_snapshot);
		const int rowSys = row(system);
		hSys += (rowSys) ? +1 : -1;
		tcs += ((rowSnap) == (rowSys)) ? +1 : -1;
		hSnap += (rowSnap) ? +1 : -1;
		tch += hSnap*hSys;
		rl.addH(hSys);
		rl.normAndSetW2(size.dimX());
		r += rl;
		// std::cerr << '\t' << (double(tcs)/size.dimX());
		m_ch += tch/size.dimX() - rl.h[0]*m_hL[yy];
		m_cs += double(tcs)/size.dimX();
	}
	// std::cerr << '\n';
	r /= size.dimY();
	m_ch /= size.dimY();
	m_cs /= size.dimY();
	return r;
}
