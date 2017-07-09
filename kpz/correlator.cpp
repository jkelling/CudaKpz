/***************************************************************************
*   Copyright 2013 - 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#include "correlator.h"

#include <cstring>
#include <vector>
#include <utility>

#include <pthread.h>

unsigned int Kpz::Correlator::s_nThreads = 0;

Kpz::Correlator::Correlator(SchedulerService* scheduler, unsigned int stop)
	: m_scheduler(scheduler),  m_snapshot(0), m_stop(stop)
{
}

Kpz::Roughness Kpz::Correlator::set(unsigned int mcs)
{
	m_mcs = mcs;
	const SystemSize& size = m_scheduler->size();
	delete m_snapshot;
	m_snapshot = new unsigned int[size.sizeW()];
	memcpy(m_snapshot, m_scheduler->system(), size.sizeW()<<2);
	int h = 0;
	Roughness r;

	switch(m_scheduler->encoding())
	{
		case SchedulerService::ENC_LOCAL:
			for (int yy=0;yy<size.dimY();++yy)
			{
				const int row = size.latticePoint(0,yy)(m_snapshot);
				h += (row & 2) ? +1 : -1;
				for (int xx=1;xx<size.dimX();++xx)
				{
					h += (size.slopeX(xx,yy)(m_snapshot)) ? +1 : -1;
					r.addH(h);
				}
				if(row & 1)
					++h;
				else --h;

				r.addH(h);
			}
			break;
		case SchedulerService::ENC_CBBIT:
			{
				const size_t quadrant = size.sizeW() >> 2;
				for (int yy=0;yy<size.dimY();++yy)
				{
					unsigned int index[2];
					index[0] = index[1] = unsigned(yy<<(size.lDimX()-6));
					index[1] += unsigned(quadrant<<1);
					const int sigy = yy&1;
					index[sigy] += quadrant; // y
					h += (Slope(index[sigy], 0)(m_snapshot) & 1) ? +1 : -1;
					index[sigy] -= quadrant; // x
					for (int xx=1;xx<size.dimX();)
					{
						const unsigned int cell[2] = {
							m_snapshot[index[sigy] + (xx>>6)], // divide xx by number of sites encoded in odd and even word combined
							m_snapshot[index[sigy^1] + (xx>>6)]
						};
						for(int x = xx&63; x < 64; ++x, ++xx)
						{
							h += ((cell[(x&1)]>>(x>>1))&1) ? +1 : -1;
							r.addH(h);
						}
					}
					if(Slope(index[sigy], 0)(m_snapshot) & 1)
						++h;
					else --h;
					r.addH(h);
				}
			}
			break;
		default:
			std::cerr << "[BUG][Correlator] Encoding notimplemented.\n";
			exit(1);
	}
	r.normAndSetW2(size.size());
	m_h = r.h[0];

	m_ch = m_cs = 1.;
	return r;
}

struct Kpz::Correlator::CorrelateData
{
	Correlator* m_correlator;
	long long m_tcs;
	double m_ch;
	Roughness m_r;
	int m_hSys, m_hSnap;
	int m_ymin, m_ysup;

	CorrelateData()
		: m_correlator(0)
		  , m_tcs(0), m_ch(0.)
	{}
	CorrelateData(Correlator* corr, int hSys = 0, int hSnap = 0)
		: m_correlator(corr)
		  , m_tcs(0), m_ch(0.)
		  , m_hSys(hSys), m_hSnap(hSnap)
		  , m_ymin(0)
	{
		m_ysup = m_correlator->m_scheduler->size().dimYW();
	}
	CorrelateData(Correlator* corr, int hSys, int hSnap, int ymin, int ysup)
		: m_correlator(corr)
		  , m_tcs(0), m_ch(0.)
		  , m_hSys(hSys), m_hSnap(hSnap)
		  , m_ymin(ymin), m_ysup(ysup)
	{}

	void correlate();
	/*! Cache optimization by processing by basic-cell */
	void correlateCO();
	void roughnessCO();
	/*! ENC_CBBIT */
	void correlateCBBit();
	/*! ENC_CBBIT */
	void roughnessCBBit();
};

void Kpz::Correlator::CorrelateData::correlate()
{
	const SystemSize& size = m_correlator->m_scheduler->size();
	const unsigned int* system = m_correlator->m_scheduler->system();
	int hSys = m_hSys, hSnap = m_hSnap;
	const int ysup = m_ysup<<L_BASIC_CELL_DIM_Y;
	for (int yy=(m_ymin<<L_BASIC_CELL_DIM_Y);yy<ysup;++yy)
	{
		const LatticePoint row = size.latticePoint(0,yy);
		const int rowSnap = row(m_correlator->m_snapshot);
		const int rowSys = row(system);
		m_tcs += ((rowSnap&2) == (rowSys&2)) ? +1 : -1;
		hSys += (rowSys&2) ? +1 : -1;
		hSnap += (rowSnap&2) ? +1 : -1;
		for (int xx=1;xx<size.dimX();++xx)
		{
			const LatticePoint col = size.latticePoint(xx,yy);
			const int colSnap = col(m_correlator->m_snapshot);
			const int colSys = col(system);
			hSys += (colSys&1) ? +1 : -1;
			m_tcs += ((colSnap&1) == (colSys&1)) ? +1 : -1;
			hSnap += (colSnap&1) ? +1 : -1;
			m_tcs += ((colSnap&2) == (colSys&2)) ? +1 : -1;
			m_r.addH(hSys);
			m_ch += hSnap*(double)hSys;
		}
		hSys += (rowSys&1) ? +1 : -1;
		m_tcs += ((rowSnap&1) == (rowSys&1)) ? +1 : -1;
		hSnap += (rowSnap&1) ? +1 : -1;
		m_ch += hSnap*(double)hSys;
		m_r.addH(hSys);
	}
}

void Kpz::Correlator::CorrelateData::correlateCO()
{
	const unsigned int* system = m_correlator->m_scheduler->system();

	int hSnapY = m_hSnap, hSysY = m_hSys;
	size_t cell = m_correlator->m_scheduler->size().indexW(0,m_ymin);
	const size_t dimXW = m_correlator->m_scheduler->size().dimXW();
	for (int yw=m_ymin;yw<m_ysup;++yw)
	{
		LatticePoint loc(0, 0);
		unsigned int locSys = system[cell];
		unsigned int locSnap = m_correlator->m_snapshot[cell];
		int hSnapX, hSysX;

#define update(hSysVar,hSnapVar) \
		m_tcs += ((snap&1) == (sys&1)) ? +1 : -1; \
		m_tcs += ((snap&2) == (sys&2)) ? +1 : -1; \
		m_r.addH(hSysVar); \
		m_ch += hSnapVar*(double)hSysVar;
#define updateX_Y(hSysVarX,hSnapVarX,hSysVarY,hSnapVarY) \
		{ \
			int sys = loc(locSys); \
			int snap = loc(locSnap); \
			hSysVarX = hSysVarY += (sys&2) ? +1 : -1; \
			hSnapVarX = hSnapVarY += (snap&2) ? +1 : -1; \
			int x = 0; \
			update(hSysVarX,hSnapVarX) \
			for(int x = 1; x < BASIC_CELL_DIM_X; ++x) { \
				loc.shift = SystemSize::shift(x,y); \
				sys = loc(locSys); \
				snap = loc(locSnap); \
				hSysVarX += (sys&1) ? +1 : -1; \
				hSnapVarX += (snap&1) ? +1 : -1; \
				update(hSysVarX,hSnapVarX) \
			} \
		}
#define updateX_X(hSysVarX,hSnapVarX,hSysVarY,hSnapVarY) \
		{ \
			int sys = loc(locSys); \
			int snap = loc(locSnap); \
			hSysVarX = hSysVarY += (sys&1) ? +1 : -1; \
			hSnapVarX = hSnapVarY += (snap&1) ? +1 : -1; \
			int x = 0; \
			update(hSysVarX,hSnapVarX) \
			for(int x = 1; x < BASIC_CELL_DIM_X; ++x) { \
				loc.shift = SystemSize::shift(x,y); \
				sys = loc(locSys); \
				snap = loc(locSnap); \
				hSysVarX += (sys&1) ? +1 : -1; \
				hSnapVarX += (snap&1) ? +1 : -1; \
				update(hSysVarX,hSnapVarX) \
			} \
		}

		int y = 0;
		updateX_Y(hSysX,hSnapX, hSysY,hSnapY); // first line, advance global X-height
		for(y = 1; y < 4; ++y)
		{
			int hSys, hSnap;
			loc.shift = SystemSize::shift(0,y);
			updateX_Y(hSys,hSnap, hSysY,hSnapY);
		}

		++cell;
		for (int xw=1;xw<dimXW;++xw, ++cell)
		{
			locSys = system[cell];
			locSnap = m_correlator->m_snapshot[cell];
			int hSnapYl = hSnapX, hSysYl = hSysX;
			y = 0;
			loc.shift = SystemSize::shift(0,0);
			updateX_X(hSysX,hSnapX, hSysYl,hSnapYl); // first line, advance global X-height
			for(y = 1; y < 4; ++y)
			{
				loc.shift = SystemSize::shift(0,y);
				int hSys, hSnap;
				updateX_Y(hSys,hSnap, hSysYl,hSnapYl);
			}
		}
	}
#undef update
#undef updateX_X
#undef updateX_Y
}

void Kpz::Correlator::CorrelateData::roughnessCO()
{
	const unsigned int* system = m_correlator->m_scheduler->system();

	int hSysY = m_hSys;
	size_t cell = m_correlator->m_scheduler->size().indexW(0,m_ymin);
	const size_t dimXW = m_correlator->m_scheduler->size().dimXW();
	for (int yw=m_ymin;yw<m_ysup;++yw)
	{
		LatticePoint loc(0, 0);
		unsigned int locSys = system[cell];
		int hSysX;

#define update(hSysVar) \
		m_r.addH(hSysVar);
#define updateX_Y(hSysVarX,hSysVarY) \
		{ \
			int sys = loc(locSys); \
			hSysVarX = hSysVarY += (sys&2) ? +1 : -1; \
			int x = 0; \
			update(hSysVarX) \
			for(int x = 1; x < BASIC_CELL_DIM_X; ++x) { \
				loc.shift = SystemSize::shift(x,y); \
				sys = loc(locSys); \
				hSysVarX += (sys&1) ? +1 : -1; \
				update(hSysVarX) \
			} \
		}
#define updateX_X(hSysVarX,hSysVarY) \
		{ \
			int sys = loc(locSys); \
			hSysVarX = hSysVarY += (sys&1) ? +1 : -1; \
			int x = 0; \
			update(hSysVarX) \
			for(int x = 1; x < BASIC_CELL_DIM_X; ++x) { \
				loc.shift = SystemSize::shift(x,y); \
				sys = loc(locSys); \
				hSysVarX += (sys&1) ? +1 : -1; \
				update(hSysVarX) \
			} \
		}

		int y = 0;
		updateX_Y(hSysX, hSysY); // first line, advance global X-height
		for(y = 1; y < 4; ++y)
		{
			int hSys;
			loc.shift = SystemSize::shift(0,y);
			updateX_Y(hSys, hSysY);
		}

		++cell;
		for (int xw=1;xw<dimXW;++xw, ++cell)
		{
			locSys = system[cell];
			int hSysYl = hSysX;
			y = 0;
			loc.shift = SystemSize::shift(0,0);
			updateX_X(hSysX, hSysYl); // first line, advance global X-height
			for(y = 1; y < 4; ++y)
			{
				loc.shift = SystemSize::shift(0,y);
				int hSys;
				updateX_Y(hSys, hSysYl);
			}
		}
	}
#undef update
#undef updateX_X
#undef updateX_Y
}

void Kpz::Correlator::CorrelateData::correlateCBBit()
{
	const unsigned int* system = m_correlator->m_scheduler->system();
	const SystemSize& size = m_correlator->m_scheduler->size();
	int hSys = m_hSys, hSnap = m_hSnap;
	const size_t quadrant = size.sizeW() >> 2;
	for (int yy=m_ymin<<2; yy<(m_ysup<<2); ++yy) // ymin and ysup are in words for ENC_LOCAL
	{
		unsigned int index[2];
		index[0] = index[1] = unsigned(yy<<(size.lDimX()-6));
		index[1] += unsigned(quadrant<<1);
		const int sigy = yy&1;
		index[sigy] += quadrant; // y
		hSys += (Slope(index[sigy], 0)(system) & 1) ? +1 : -1;
		hSnap += (Slope(index[sigy], 0)(m_correlator->m_snapshot) & 1) ? +1 : -1;
		index[sigy] -= quadrant; // x
		for (int xx=1;xx<size.dimX();)
		{
			const unsigned int cellSys[2] = {
				// divide xx by number of sites encoded in odd and even word combined
				system[index[sigy] + (xx>>6)],
				system[index[sigy^1] + (xx>>6)]
			};
			const unsigned int cellSysY[2] = {
				system[index[sigy] + (xx>>6) + quadrant],
				system[index[sigy^1] + (xx>>6) + quadrant]
			};
			const unsigned int cellSnap[2] = {
				m_correlator->m_snapshot[index[sigy] + (xx>>6)],
				m_correlator->m_snapshot[index[sigy^1] + (xx>>6)]
			};
			const unsigned int cellSnapY[2] = {
				m_correlator->m_snapshot[index[sigy] + (xx>>6) + quadrant],
				m_correlator->m_snapshot[index[sigy^1] + (xx>>6) + quadrant]
			};
			for(int x = xx&63; x < 64; ++x, ++xx)
			{
				{
					const int sXSys = ((cellSys[(x&1)]>>(x>>1))&1);
					const int sXSnap = ((cellSnap[(x&1)]>>(x>>1))&1);
					hSys += sXSys ? +1 : -1;
					hSnap += sXSnap ? +1 : -1;
					m_tcs += (sXSys == sXSnap) ? +1 : -1;
				}
				m_tcs += (((cellSysY[(x&1)]>>(x>>1))&1) == ((cellSnapY[(x&1)]>>(x>>1))&1)) ? +1 : -1;
				m_r.addH(hSys);
				m_ch += hSnap*(double)hSys;
			}
		}
		{
			const int sXSys = (Slope(index[sigy], 0)(system) & 1);
			const int sXSnap = (Slope(index[sigy], 0)(m_correlator->m_snapshot) & 1);
			hSys += sXSys ? +1 : -1;
			hSnap += sXSnap ? +1 : -1;
			m_tcs += (sXSys == sXSnap) ? +1 : -1;
		}
		// corr of y slopes
		m_tcs += ((Slope(index[sigy]+quadrant, 0)(system) & 1)
				== (Slope(index[sigy] + quadrant, 0)(m_correlator->m_snapshot) & 1)) ? +1 : -1;
		m_r.addH(hSys);
		m_ch += hSnap*(double)hSys;
	}
}

void Kpz::Correlator::CorrelateData::roughnessCBBit()
{
	const unsigned int* system = m_correlator->m_scheduler->system();
	const SystemSize& size = m_correlator->m_scheduler->size();
	int h = m_hSys;
	const size_t quadrant = size.sizeW() >> 2;
	for (int yy=m_ymin<<2; yy<(m_ysup<<2); ++yy) // ymin and ysup are in words for ENC_LOCAL
	{
		unsigned int index[2];
		index[0] = index[1] = unsigned(yy<<(size.lDimX()-6));
		index[1] += unsigned(quadrant<<1);
		const int sigy = yy&1;
		index[sigy] += quadrant; // y
		h += (Slope(index[sigy], 0)(system) & 1) ? +1 : -1;
		index[sigy] -= quadrant; // x
		for (int xx=1;xx<size.dimX();)
		{
			const unsigned int cell[2] = {
				system[index[sigy] + (xx>>6)], // divide xx by number of sites encoded in odd and even word combined
				system[index[sigy^1] + (xx>>6)]
			};
			for(int x = xx&63; x < 64; ++x, ++xx)
			{
				h += ((cell[(x&1)]>>(x>>1))&1) ? +1 : -1;
				m_r.addH(h);
			}
		}
		if(Slope(index[sigy], 0)(system) & 1)
			++h;
		else --h;
		m_r.addH(h);
	}
}

Kpz::Roughness Kpz::Correlator::correlate()
{
	CorrelateData corr(this);
	switch(m_scheduler->encoding())
	{
		case SchedulerService::ENC_LOCAL:
			corr.correlateCO();
			break;
		case SchedulerService::ENC_CBBIT:
			corr.correlateCBBit();
			break;
		default:
			std::cerr << "[BUG][Correlator] Encoding notimplemented.\n";
			exit(1);
	}

	const SystemSize& size = m_scheduler->size();
	corr.m_r.normAndSetW2(size.size());
	m_ch = corr.m_ch/size.size() - corr.m_r.h[0]*m_h;
	m_cs = double(corr.m_tcs)/size.size()/2.;
	return corr.m_r;
}

Kpz::Roughness Kpz::Correlator::roughness()
{
	CorrelateData corr(this);
	switch(m_scheduler->encoding())
	{
		case SchedulerService::ENC_LOCAL:
			corr.roughnessCO();
			break;
		case SchedulerService::ENC_CBBIT:
			corr.roughnessCBBit();
			break;
		default:
			std::cerr << "[BUG][Correlator] Encoding notimplemented.\n";
			exit(1);
	}

	const SystemSize& size = m_scheduler->size();
	corr.m_r.normAndSetW2(size.size());
	return corr.m_r;
}

void* Kpz::Correlator::correlateThreadCO(void* cd)
{
	auto data = reinterpret_cast<Correlator::CorrelateData*>(cd);
	data->correlateCO();
}

void* Kpz::Correlator::correlateThreadCBBit(void* cd)
{
	auto data = reinterpret_cast<Correlator::CorrelateData*>(cd);
	data->correlateCBBit();
}

Kpz::Roughness Kpz::Correlator::correlateThreaded()
{
	if(!s_nThreads)
		return correlate();

	void* (*corrFct)(void*);
	switch(m_scheduler->encoding())
	{
		case SchedulerService::ENC_LOCAL:
			corrFct = correlateThreadCO;
			break;
		case SchedulerService::ENC_CBBIT:
			corrFct = correlateThreadCBBit;
			break;
		default:
			std::cerr << "[BUG][Correlator] Encoding notimplemented.\n";
			exit(1);
	}

	std::vector<std::pair<pthread_t, CorrelateData> > threads(s_nThreads);
	const unsigned int* system = m_scheduler->system();
	const SystemSize& size = m_scheduler->size();
	const size_t div = size.dimYW()/s_nThreads;
	const size_t mod = size.dimYW()%s_nThreads;

	int hSys = 0, hSnap = 0;
	int ymin = 0, ysup = div + (mod>0);
	size_t cell = 0;
	unsigned int t = 0;

	threads[t].second.m_ymin = ymin;
	threads[t].second.m_ysup = ysup;
	threads[t].second.m_hSys = hSys;
	threads[t].second.m_hSnap = hSnap;
	threads[t].second.m_correlator = this;
	pthread_create(&threads[t].first, NULL, corrFct, &threads[t].second);

	for(t=1; t < s_nThreads; ++t)
	{
		switch(m_scheduler->encoding())
		{
			case SchedulerService::ENC_LOCAL:
				for(; ymin < ysup; ++ymin, cell+=size.dimXW())
				{
					const unsigned int locSys = system[cell];
					const unsigned int locSnap = m_snapshot[cell];
					LatticePoint loc(0,0);
					hSys += (loc(locSys)&2) ? +1 : -1;
					hSnap += (loc(locSnap)&2) ? +1 : -1;
					for(int y = 1; y < 4; ++y)
					{
						loc.shift = SystemSize::shift(0,y);
						hSys += (loc(locSys)&2) ? +1 : -1;
						hSnap += (loc(locSnap)&2) ? +1 : -1;
					}
				}
				break;
			case SchedulerService::ENC_CBBIT:
				{
					int y = ymin << 2;
					unsigned int index[2];
					index[0] = index[1] = unsigned(y<<(size.lDimX()-6))+(size.sizeW()>>2);
					index[1] += size.sizeW()>>1;
					index[(y&1)^1] += 1<<(size.lDimX()-6);
					for(; y < ysup<<2; ++y)
					{
						hSys += (system[index[y&1]] & 1) ? +1 : -1;
						hSnap += (m_snapshot[index[y&1]] & 1) ? +1 : -1;
						index[y&1] += 1<<(size.lDimX()-5);
					}
					ymin = ysup;
				}
				break;
		}

		threads[t].second.m_ymin = ymin;
		ysup += div + (mod > t);
		threads[t].second.m_ysup = ysup;
		threads[t].second.m_hSys = hSys;
		threads[t].second.m_hSnap = hSnap;
		threads[t].second.m_correlator = this;
		pthread_create(&threads[t].first, NULL, corrFct, &threads[t].second);
	}

	t = 0;
	pthread_join(threads[t].first, NULL);
	Roughness r(threads[t].second.m_r);
	m_ch = threads[t].second.m_ch;
	long long tcs = threads[t].second.m_tcs;

	for(t = 1; t < s_nThreads; ++t)
	{
		pthread_join(threads[t].first, NULL);
		r += threads[t].second.m_r;
		m_ch += threads[t].second.m_ch;
		tcs += threads[t].second.m_tcs;
	}

	r.normAndSetW2(size.size());
	m_ch = m_ch/size.size() - r.h[0]*m_h;
	m_cs = double(tcs)/size.size()/2.;
	return r;
}

void* Kpz::Correlator::roughnessThreadCO(void* cd)
{
	auto data = reinterpret_cast<Correlator::CorrelateData*>(cd);
	data->roughnessCO();
}

void* Kpz::Correlator::roughnessThreadCBBit(void* cd)
{
	auto data = reinterpret_cast<Correlator::CorrelateData*>(cd);
	data->roughnessCBBit();
}

Kpz::Roughness Kpz::Correlator::roughnessThreaded()
{
	if(!s_nThreads)
		return roughness();

	void* (*roughnessFct)(void*);
	switch(m_scheduler->encoding())
	{
		case SchedulerService::ENC_LOCAL:
			roughnessFct = &roughnessThreadCO;
			break;
		case SchedulerService::ENC_CBBIT:
			roughnessFct = &roughnessThreadCBBit;
			break;
		default:
			std::cerr << "[BUG][Correlator] Encoding notimplemented.\n";
			exit(1);
	}

	std::vector<std::pair<pthread_t, CorrelateData> > threads(s_nThreads);
	const unsigned int* system = m_scheduler->system();
	const SystemSize& size = m_scheduler->size();
	const size_t div = size.dimYW()/s_nThreads;
	const size_t mod = size.dimYW()%s_nThreads;

	int hSys = 0;
	int ymin = 0, ysup = div + (mod>0);
	size_t cell = 0;
	unsigned int t = 0;

	threads[t].second.m_ymin = ymin;
	threads[t].second.m_ysup = ysup;
	threads[t].second.m_hSys = hSys;
	threads[t].second.m_correlator = this;
	pthread_create(&threads[t].first, NULL, roughnessFct, &threads[t].second);

	for(t=1; t < s_nThreads; ++t)
	{
		switch(m_scheduler->encoding())
		{
			case SchedulerService::ENC_LOCAL:
				for(; ymin < ysup; ++ymin, cell+=size.dimXW())
				{
					const unsigned int locSys = system[cell];
					LatticePoint loc(0,0);
					hSys += (loc(locSys)&2) ? +1 : -1;
					for(int y = 1; y < 4; ++y)
					{
						loc.shift = SystemSize::shift(0,y);
						hSys += (loc(locSys)&2) ? +1 : -1;
					}
				}
				break;
			case SchedulerService::ENC_CBBIT:
				{
					int y = ymin << 2;
					unsigned int index[2];
					index[0] = index[1] = unsigned(y<<(size.lDimX()-6))+(size.sizeW()>>2);
					index[1] += size.sizeW()>>1;
					index[(y&1)^1] += 1<<(size.lDimX()-6);
					for(; y < ysup<<2; ++y)
					{
						hSys += (system[index[y&1]] & 1) ? +1 : -1;
						index[y&1] += 1<<(size.lDimX()-5);
					}
					ymin = ysup;
				}
				break;
		}

		threads[t].second.m_ymin = ymin;
		ysup += div + (mod > t);
		threads[t].second.m_ysup = ysup;
		threads[t].second.m_hSys = hSys;
		threads[t].second.m_correlator = this;
		pthread_create(&threads[t].first, NULL, roughnessFct, &threads[t].second);
	}

	t = 0;
	pthread_join(threads[t].first, NULL);
	Roughness r(threads[t].second.m_r);

	for(t = 1; t < s_nThreads; ++t)
	{
		pthread_join(threads[t].first, NULL);
		r += threads[t].second.m_r;
	}

	r.normAndSetW2(size.size());
	return r;
}

void Kpz::Correlator::setNThreads(unsigned int n)
{
	if(n > 1)
		s_nThreads = n;
	else
		s_nThreads = 0;
}

#include <KMCsplash.h>

void Kpz::Correlator::save(splash::DataCollector* data, const std::string& prefix)
{
	const std::string sysname = prefix + ".system";
	splash::Selection sel(m_scheduler->size().splashDimensions());
	data->write(m_mcs, KmcSplash::ColTypeUInt32, 2, sel, sysname.c_str(), m_snapshot);
	data->writeAttribute(m_mcs, KmcSplash::ColTypeUInt32, 0, (prefix + ".corrTag").c_str(), &m_stop);
	data->writeAttribute(m_mcs, KmcSplash::ColTypeDouble, sysname.c_str(), "avgH", &m_h);
}

Kpz::Correlator* Kpz::Correlator::fromH5(splash::DataCollector* data, int id, SchedulerService* scheduler, const std::string& prefix)
{
	const std::string sysname = prefix + ".system";

	Correlator* corr = 0;
	int stop;
	try {
		data->readAttribute(id, 0, (prefix + ".corrTag").c_str(), &stop);
	}
	catch(splash::DCException e)
	{
		// try legacy libSplash/hdf5
		try {
			data->readAttribute(id, 0, (prefix + "corrTag").c_str(), &stop);
		}
		catch(splash::DCException e)
		{
			return 0;
		}
	}
	corr = new Correlator(scheduler, stop);
	try {
		data->readAttribute(id, sysname.c_str(), "avgH", &corr->m_h);
		corr->m_mcs = id;
		corr->m_snapshot = scheduler->size().init();
		auto dim = scheduler->size().splashDimensions();
		data->read(id, sysname.c_str(), dim, corr->m_snapshot);
		return corr;
	}
	catch(splash::DCException e)
	{
		std::cerr << "[EE][Correlator::fromH5] " << e.what() << '\n';
		return 0;
	}
}

std::ostream& operator<< (std::ostream& o, const Kpz::Correlator& c) {
	return o << c.mcs() << '\t' << c.ch() << '\t' << c.cs();
}
