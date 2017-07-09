/***************************************************************************
*   Copyright 2011 - 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#include "schedulerService.h"

#include "kpzConst.h"
#include "lattice.h"
#include "correlator.h"
#include "correlatorLines.h"

#include <benchmarking.h>
#include <kmcMath.h>

#include <cstring>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <sstream>
//#include <cmath>
#include <iomanip>
#include <utility>
#include <functional>
#include <stdexcept>

#include <mpi.h>

using namespace Kmc;

bool Kpz::SchedulerService::initDisorder2()
{
	if(!m_disorder)
		m_disorder = m_size.init();
	return m_disorder != 0;
}
bool Kpz::SchedulerService::generateDisorder2(double fill)
{
	if(!initDisorder2())
		return false;
	for(int a = 0; a < m_size.sizeW(); ++a)
	{
		int t = 0;
		for(int s = 0; s < 32; s+=2)
		{
			if(dsfmt_genrand_close_open(&m_dsfmt) < fill)
			{
				t |= (1<<s);
			}
		}
		m_disorder[a] = t;
	}
	return true;
}

std::ostream& operator<<(std::ostream& o, const Kpz::Roughness& r)
{
	o << r.w2;
	for(auto a : r.h)
		o << '\t' << a;
	return o;
}

bool Kpz::SchedulerService::loadXYZ(const char* file)
{
	return true;
}

bool Kpz::SchedulerService::loadBit(const char* file)
{
	return true;
}

Kpz::SchedulerService::SchedulerService(int lx, int ly, Encoding enc)
	: m_system(0), m_size(lx, ly), out(&std::cout), m_roughnessThread(0)
	, m_mcsThread(0)
	, m_mcs(0), m_disorder(0), m_deviceMcs(0)
	, m_silentRoughness(false), m_outPrecision(12)
	, m_callMcs(&callMcs)
	, m_encoding(enc)
	, m_fAnalysisOnGPU(false)
{
	initURandom();

	if(!setSize())
		exit(1);
}

Kpz::SchedulerService::SchedulerService(const SystemSize& size, Encoding enc)
	: m_system(0), m_size(size), m_roughnessThread(0), out(&std::cout)
	, m_mcsThread(0)
	, m_mcs(0), m_disorder(0), m_deviceMcs(0)
	, m_silentRoughness(false), m_outPrecision(12)
	, m_callMcs(&callMcs)
	, m_encoding(enc)
	, m_fAnalysisOnGPU(false)
{
	initURandom();

	if(!setSize())
		exit(1);
}

Kpz::SchedulerService::SchedulerService(SchedulerService &&other)
	: m_system(other.m_system), m_size(other.m_size), m_roughnessThread(0), out(other.out)
	, m_mcsThread(0)
	, m_mcs(other.m_mcs), m_deviceMcs(other.m_deviceMcs)
	, m_disorder(other.m_disorder)
	, m_silentRoughness(other.m_silentRoughness), m_outPrecision(other.m_outPrecision)
	, m_dsfmt(other.m_dsfmt)
	, m_callMcs(other.m_callMcs)
	, m_encoding(other.encoding())
	, m_fAnalysisOnGPU(false)
{
	other.m_system = 0;
	other.m_disorder = 0;
}

Kpz::SchedulerService::~SchedulerService()
{
	cout().flush();
	delete[] m_system;
	delete[] m_disorder;
	for(typename std::list<Correlator*>::iterator a = m_correlateTags.begin();
			a != m_correlateTags.end(); a = m_correlateTags.erase(a))
		delete *a;

	std::cout << "# When publishing data generated using CudaKPZ, or referring to it, plese cite:\n"
		"#\tKelling, Jeffrey, Ódor, Géza, Gemming, Sibylle: Suppressing correlations in massively"
		"parallel simulations of lattice models, https://arxiv.org/abs/1705.01022, 2017\n"
		"#\tKelling, Jeffrey, Ódor, Géza, Gemming, Sibylle: Bit-Vectorized GPU"
		"Implementation of a Stochastic Cellular Automaton Model for Surface Growth,"
		"2016 IEEE International Conference on Intelligent Engineering Systems, 2016."
		"INES '16, IEEE, https://doi.org/10.1109/INES.2016.7555127, 2016\n"
		"# The former describes implementations of RS dynamics and the latter those for SCA dynamics."
		<< std::endl;
}

void Kpz::SchedulerService::initHomogeneous(Encoding encoding)
{
	switch(encoding)
	{
		case ENC_LOCAL:
			for(int a = 0; a < m_size.sizeW(); ++a)
				m_system[a] = 0x33CC33CC;
			break;
		case ENC_CBBIT:
			// all sublattice 0 x = 0 and 0 y = 0
			memset(m_system, 0, (m_size.sizeW() << 1));
			// all sublattice 1 x = 1 and 1 y = 1
			memset(m_system + (m_size.sizeW() >> 1), 0xFF, (m_size.sizeW() << 1));
			break;
		default:
			std::cerr << "[SchedulerService][initHomogeneous][BUG] Encoding not implemented.\n";
			exit(1);
	}
	m_encoding = encoding;
}

void Kpz::SchedulerService::nullSystem()
{
	memset(m_system, 0, m_size.sizeW()<<2);
}

void* Kpz::SchedulerService::callRoughness(void* scheduler)
{
	timeAinit
	timeAstart
	using ::operator<<;
	SchedulerService* s = reinterpret_cast<SchedulerService* >(scheduler);
	s->m_roughness = s->roughness();
	if(!s->m_silentRoughness)
	{
		std::ostringstream text;
		text << std::setprecision(s->m_outPrecision) << s->m_mcs << '\t' << s->m_roughness << '\n';
		s->cout() << text.str();
		s->cout().flush();
	}
	timeAstop("calcRoughnessAsync");
	return 0;
}

void* Kpz::SchedulerService::threadCorrelateTag(void* scheduler)
{
	using ::operator<<;
	SchedulerService* s = reinterpret_cast<SchedulerService* >(scheduler);
	s->m_correlateTags.push_back(new Correlator(s, s->m_stopCorr));
	s->m_roughness = s->m_correlateTags.back()->set(s->m_mcs);
	if(!s->m_silentRoughness)
	{
		std::ostringstream text;
		text << std::setprecision(s->m_outPrecision) << s->m_mcs << '\t' << s->m_roughness << "\tset\tcorr\ttag\n";
		s->cout() << text.str();
		s->cout().flush();
	}
	return 0;
}

void* Kpz::SchedulerService::threadCorrelate(void* scheduler)
{
	timeAinit;
	timeAstart;
	using ::operator<<;
	SchedulerService* s = reinterpret_cast<SchedulerService* >(scheduler);
	std::ostringstream text;
	bool corr = false;
	for(typename std::list<Correlator*>::iterator a = s->m_correlateTags.begin(); a != s->m_correlateTags.end();)
	{
		// std::cerr << s->m_mcs;
		if(s->m_mcs > (*a)->stop())
			a = s->m_correlateTags.erase(a);
		else
		{
			s->m_roughness = (*a)->correlateThreaded();
			if(!s->m_silentRoughness)
				text << std::setprecision(s->m_outPrecision) << s->m_mcs << '\t' << s->m_roughness << '\t' << (**a) << '\n';
			++a;
			corr = true;
		}
	}
	if(!corr)
	{
		s->m_roughness = s->roughness();
		if(!s->m_silentRoughness)
			text << std::setprecision(s->m_outPrecision) << s->m_mcs << '\t' << s->m_roughness << '\n';
	}
	if(!s->m_silentRoughness)
	{
		s->cout() << text.str();
		s->cout().flush();
	}
	timeAstop("SchedulerService::threadCorrelate");
	return 0;
}

void* Kpz::SchedulerService::threadCorrelateLinesTag(void* scheduler)
{
	using ::operator<<;
	SchedulerService* s = reinterpret_cast<SchedulerService* >(scheduler);
	s->m_correlateTags.push_back(new CorrelatorLines(s));
	std::ostringstream text;
	const Kpz::Roughness r = s->m_correlateTags.back()->set(s->m_mcs);
	text << std::setprecision(s->m_outPrecision) << s->m_mcs << '\t' << r << "\tset\tcorr\ttag\n";
	s->cout() << text.str();
	s->cout().flush();
	return 0;
}

void Kpz::SchedulerService::roughnessAsync()
{
	joinRoughness();
	pthread_create(&m_roughnessThread, NULL, &callRoughness, (void*) this);
#if 0
	timeAinit
	timeAstart
	using ::operator<<;
	cout() << mcs << '\t' << roughness() << '\n';
	timeAstop("calcRoughness");
#endif
}

void Kpz::SchedulerService::correlateTag(unsigned int stop)
{
	joinRoughness();
	m_stopCorr = stop;
	pthread_create(&m_roughnessThread, NULL, &threadCorrelateTag, (void*) this);
}

void Kpz::SchedulerService::correlate()
{
	joinRoughness();
	pthread_create(&m_roughnessThread, NULL, &threadCorrelate, (void*) this);
}

Kpz::Roughness Kpz::SchedulerService::roughness()
{
#if 1
	Correlator tmp(this);
	return tmp.roughnessThreaded();
#else
	int h = 0;
	Roughness r;
	for (int yy=0;yy<m_size.dimY();++yy)
	{
		h += (m_size.slopeY(0,yy)(m_system))?(+1):(-1);
		for (int xx=1;xx<m_size.dimX();++xx) 
		{
			h += (m_size.slopeX(xx,yy)(m_system))?(+1):(-1);
			r.h += h;
			r.h2 += h*h;
		}
		h += (m_size.slopeX(0,yy)(m_system))?(+1):(-1);
		r.h += h;
		r.h2 += h*h;
	}
	r.h /= m_size.size();
	r.h2 /= m_size.size();
	r.w2 = r.h2 - r.h*r.h;
	return r;
#endif
}

Kpz::Roughness Kpz::SchedulerService::roughnessLines()
{
	Roughness r;
	for (int yy=0;yy<m_size.dimY();++yy)
	{
		long double h = 0;
		double rh[Roughness::H_ORDERS];
		for(auto& a : rh)
			a = 0;
		for (int xx=1;xx<m_size.dimX();++xx) 
		{
			h += (m_size.slopeX(xx,yy)(m_system))?(+1):(-1);
			auto acc = 1.;
			for(auto& a : rh)
				a += (acc*=h);
		}
		h += (m_size.slopeX(0,yy)(m_system))?(+1):(-1);
		long double acc = 1.;
		for(auto& a : rh)
			a += (acc*=h);
		for(int i = 0; i < Roughness::H_ORDERS; ++i)
			r.h[i] += (rh[i] /= m_size.dimX());
		r.w2 += rh[0] - rh[0]*rh[0];
	}
	for(auto& a : r.h)
		a /= m_size.dimY();
	r.w2 /= m_size.dimY();
	return r;
}

bool Kpz::SchedulerService::setSize(int lx, int ly)
{
	if((lx == m_size.lDimX()) && (ly == m_size.lDimY()))
		return true;
	m_size.set(lx, ly);
	return setSize();
}

bool Kpz::SchedulerService::setSize()
{
	delete m_disorder; // disorder invalid
	m_disorder = 0;
	delete m_system;
	m_system = 0;
	m_system = m_size.init();
	if(!m_system)
	{
		std::cerr << "Insufficient hostmemory to create sytem of size " << m_size.lDimX() << ',' << m_size.lDimY() << " .\n";
		return false;
	}

	changedSizeEvent(); //give derived class opportunity to adjust

	std::cout << "\nSet system size to " << m_size.lDimX() << ',' << m_size.lDimY() << '\n' << m_size;
	return true;
}

void Kpz::SchedulerService::pushSystem()
{
	join();
	m_deviceMcs = m_mcs;
}

void Kpz::SchedulerService::popSystem()
{
	join();
	m_mcs = m_deviceMcs;
}

bool Kpz::SchedulerService::saveBit(const char* file, const char* meta)
{
	std::ofstream out(file);
	out << m_size.lDimX() << '\n' << m_size.lDimY() << '\n';
	if(meta)
		out << meta;
	out << '\0';
	out.close();
	out.open(file, std::ios_base::binary | std::ios_base::app);
	out.write((char*)m_system, m_size.sizeW()*sizeof(int));
	out.close();
	return true;
}

template<class SlopeX, class SlopeY>
void printXYZHeightmap(std::ostream& o, Kpz::SchedulerService* pthis, SlopeX slopeX, SlopeY slopeY)
{
	const char ELEM[] = "Fe";

	double h = 0;
	for (int yy=0;yy<pthis->size().dimY();++yy)
	{
		h += (slopeY(0,yy)(pthis->system()))?(+1):(-1);
		for (int xx=1;xx<pthis->size().dimX();++xx)
		{
			h += (slopeX(xx,yy)(pthis->system()))?(+1):(-1);
			o << ELEM << '\t' << (xx) << '\t' << (yy) << '\t' << h << '\t' << h << '\n';
		}
		h += (slopeX(0,yy)(pthis->system()))?(+1):(-1);
		o << ELEM << '\t' << (0) << '\t' << (yy) << '\t' << h << '\t' << h << '\n';
	}
}

bool Kpz::SchedulerService::saveXYZ(const char* filename)
{
	int atoms = 4 + m_size.size();

	std::ofstream o (filename);
	o << atoms << "\n\n";

	const int dimZ = std::max(m_size.dimX(), m_size.dimY());
	o << "H\t0\t0\t0\t0\n"
		<< "H\t" << (m_size.dimX()) << '\t' << 0 << '\t' << 0 << "\t0\n"
		<< "H\t" << (m_size.dimX()) << '\t' << (m_size.dimY()) << '\t' << 0 << "\t0\n"
		<< "H\t" << 0 << '\t' << (m_size.dimY()) << '\t' << 0 << "\t0\n";
	
	using namespace std::placeholders;
	switch(m_encoding)
	{
		case ENC_LOCAL:
			printXYZHeightmap(o, this
					, std::bind(&SystemSize::slopeX, &m_size, _1, _2), std::bind(&SystemSize::slopeY, &m_size, _1, _2));
			break;
		case ENC_CBBIT:
			printXYZHeightmap(o, this
					, std::bind(&SystemSize::slopeX_CBBit, &m_size, _1, _2), std::bind(&SystemSize::slopeY_CBBit, &m_size, _1, _2));
			break;
		default:
			std::cerr << "[SchedulerService][saveXYZ][BUG] Encoding not implemented.\n";
			exit(1);
	}

	o.close();
	return true;
}

bool Kpz::SchedulerService::load(const char* file)
{
	return true;
}

void Kpz::SchedulerService::reSeed(int seed)
{
	if(seed)
		dsfmt_init_gen_rand(&m_dsfmt, seed);
	else
		dsfmt_init_gen_rand(&m_dsfmt, rand());
	for( int a = 0; a < 100000; ++a) // arbitray large number 
	{
		dsfmt_genrand_close_open(&m_dsfmt);
	}
}

void Kpz::SchedulerService::randomize()
{

}

void Kpz::SchedulerService::supplySeed(int std, int dsfmt)
{
	srand(std);
	reSeed(dsfmt);
	std::cout << "[SchedulerService][supplySeed] srand= " << std << " seed= " << dsfmt << '\n';
}

void Kpz::SchedulerService::initURandom()
{
	//initialize RNG using /dev/urandom
	std::ifstream r("/dev/urandom", std::ios_base::binary);
	int random[2];
	r.get((char*)random, 8);
	r.close();
	srand(random[0]^=time(0));

	reSeed(random[1]^=time(0)^rand());
	std::cout << "[SchedulerService][initURandom] srand= " << random[0] << " seed= " << random[1] << '\n';
}

bool Kpz::SchedulerService::copySystem(const SchedulerService* s)
{
	m_size = s->size();
	if(!setSize())
		return false;
	memcpy(m_system, s->system(), m_size.sizeW()<<2);
	return true;
}

bool Kpz::SchedulerService::copyDisorder(const SchedulerService* s)
{
	if(m_size.lDimX() == s->size().lDimX() && m_size.lDimY() == s->size().lDimY())
	{
		if(!m_disorder)
			if(!initDisorder2())
				return false;
		memcpy(m_disorder, s->disorder(), m_size.sizeW()<<2);
		return true;
	}
	return false;
}

Kpz::SchedulerService::RngPackage::RngPackage(dsfmt_t& dsfmt)
{
	for(int a = 0; a < RNG_SYNC_PACAGE; ++a)
		m[a] = dsfmt_genrand_close_open(&dsfmt);
}

#include <KMCsplash.h>

splash::DataCollector* Kpz::SchedulerService::writeH5(const char* file, const std::string& prefix)
{
	auto data = new splash::SerialDataCollector(1);
	splash::DataCollector::FileCreationAttr fcattr;
	splash::DataCollector::initFileCreationAttr(fcattr);
	fcattr.enableCompression = true;
	if(MPI::Is_initialized())
	{
		fcattr.mpiSize = splash::Dimensions(MPI::COMM_WORLD.Get_size(), 1, 1);
		fcattr.mpiPosition = splash::Dimensions(MPI::COMM_WORLD.Get_rank(), 0, 0);
	}
	data->open(file, fcattr);

	try {
	if(!collectH5(data, prefix))
	{
		closeH5(data);
		return 0;
	}
	} catch (splash::DCException e) {
		std::cerr << "[SchedulerService][writeH5][DBG] closeH5(), before forwarding exception.\n";
		closeH5(data);
		throw e;
	}

	return data;
}

splash::DataCollector* Kpz::SchedulerService::readH5(const char* file, const std::string& prefix)
{
	auto data = readH5open(file);

	if(!retrieveH5(data, prefix))
	{
		closeH5(data);
		return 0;
	}

	return data;
}

splash::DataCollector* Kpz::SchedulerService::readH5open(const char* file)
{
	auto data = new splash::SerialDataCollector(1);
	splash::DataCollector::FileCreationAttr fcattr;
	splash::DataCollector::initFileCreationAttr(fcattr);
	if(MPI::Is_initialized())
	{
		fcattr.mpiSize = splash::Dimensions(MPI::COMM_WORLD.Get_size(), 1, 1);
		fcattr.mpiPosition = splash::Dimensions(MPI::COMM_WORLD.Get_rank(), 0, 0);
	}
	fcattr.fileAccType = splash::DataCollector::FAT_READ;
	data->open(file, fcattr);
	return data;
}

bool Kpz::SchedulerService::collectH5(splash::DataCollector* data, const std::string& prefix)
{
	splash::ColTypeDouble4Array ctPQ;
	data->writeGlobalAttribute(ctPQ, (prefix + ".disorderP").c_str(), m_disorderP);
	data->writeGlobalAttribute(ctPQ, (prefix + ".disorderQ").c_str(), m_disorderQ);
	auto s = m_size.splashDimensions();
	splash::Selection sel(m_size.splashDimensions());

	// std::ostringstream tmpout;
	// auto printDim = [&tmpout](const splash::Dimensions& dim, const char* name) -> std::ostream& {
	// 	return std::cerr << "[SchedulerService][collectH5][DBG]\t" << name << ": "
	// 		<< dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
	// };
	// tmpout << "[SchedulerService][collectH5][DBG] splash::Selection:\n";
	// printDim(sel.size, "size");
	// printDim(sel.count, "count");
	// printDim(sel.offset, "offset");
	// printDim(sel.stride, "stride");
	// std::cerr << tmpout.str();

	data->write(m_mcs, KmcSplash::ColTypeUInt32, 2, sel, (prefix + ".system").c_str(), m_system);
	if(m_disorder)
		data->write(m_mcs, KmcSplash::ColTypeUInt32, 2, sel, (prefix + ".disorder").c_str(), m_disorder);

	for(auto a = m_correlateTags.begin(); a != m_correlateTags.end(); ++a)
		(*a)->save(data, prefix);

	KmcSplash::writeDSFMT(data, m_mcs, &m_dsfmt, prefix);

	return true;
}

bool Kpz::SchedulerService::retrieveH5(splash::DataCollector* data, const std::string& prefix)
{
	size_t nids = 0;
	data->getEntryIDs(0, &nids);
	int ids[nids];
	data->getEntryIDs(ids, &nids);
	std::sort(ids,&ids[nids]);

	m_mcs = ids[nids-1];
	std::cout << "[readH5] Loading Hdf5 file at t = " << m_mcs << " MCS\n";

	splash::Dimensions size;
	std::string sysname = prefix + ".system";
	data->read(m_mcs, sysname.c_str(), size, 0);
	int lx = log2Sec(size[0]);
	int ly = log2Sec(size[1]);
	if(lx < 2 || ly < 2)
	{
		std::cerr << "[readH5][EE] Hdf5 file contains invalid system size " << size[0] << ' ' << size[1] << " .\n";
		return false;
	}
	setSize(lx+2,ly+2);
	data->read(m_mcs, sysname.c_str(), size, m_system);

	sysname = prefix + ".disorder";
	try {
		data->read(m_mcs, sysname.c_str(), size, 0);
		if(lx != log2Sec(size[0]) || ly != log2Sec(size[1]))
		{
			std::cerr << "[readH5][WW] Hdf5 file contains disorder. but of incompatible size " << size[0] << ' ' << size[1]
				<< " . Going to ignore it.\n";
		}
		else
		{
			if(!initDisorder2())
			{
				std::cerr << "[readH5][EE] insufficent host memory to read disorder.\n";
				return false;
			}
			data->read(m_mcs, sysname.c_str(), size, m_disorder);
		}
		std::cout << "[readH5] Read disorder from Hdf5 file.\n";
	}
	catch(splash::DCException) {
	}
	//\todo what about disorderP and disorderQ ?

	for(int a = nids-2; a >= 0; --a)
	{
		auto corr = Correlator::fromH5(data, ids[a], this, prefix);
		if(corr)
		{
			m_correlateTags.push_back(corr);
			std::cout << "[readH5] found corrTag at " << corr->mcs() << ", stop at " << corr->stop() << '\n';
		}
	}

	if(KmcSplash::readDSFMT(data, m_mcs, &m_dsfmt, prefix))
		std::cout << "[readH5] restored dSFMT state from file.\n";

	return true;
}

int Kpz::SchedulerService::mcsH5(splash::DataCollector* data)
{
	size_t nids = 0;
	data->getEntryIDs(0, &nids);
	int ids[nids];
	data->getEntryIDs(ids, &nids);
	return *std::max_element(ids,&ids[nids]);
}

void Kpz::SchedulerService::closeH5(splash::DataCollector* data)
{
	data->close();
	delete data;
}

// async mcs* calls
void* Kpz::SchedulerService::callMcs(void* vargs)
{
	timeAinit;
	timeAstart;
	auto args = reinterpret_cast<CallMcsArgs*>(vargs);
	(args->scheduler->*(args->fkt))(args->n);
	timeAstopS("SchedulerService::callMcs_" << args->n);
	return 0;
}

int Kpz::SchedulerService::syncRNGMPI(MPI::Comm& comm)
{
	static const int root = 0;
	std::ostringstream ss;
	// ss << "Rank " << comm.Get_rank() << '\n';
	// ss << "\trand before sync " << genrand_close_open() << '\n';
	comm.Bcast(&m_dsfmt, sizeof(dsfmt_t), MPI_CHAR, root);
	// ss << "\trand after sync " << genrand_close_open() << '\n';
	// std::cerr << ss.str();
	return 0;
}

void Kpz::SchedulerService::setMcs(unsigned int mcs)
{
	m_mcs = mcs;
}

bool Kpz::SchedulerService::consitencyCheckThread(std::pair<std::ostream&, std::mutex>* out, const unsigned int* system, std::vector<int>& sumY
		, int minTidx, int supTidx, int ydivW, int ymodW) const
{
	bool ret = true;
	if(minTidx < supTidx-1)
	{
		int workLeft = (supTidx-minTidx)/2;
		std::vector<int> leftSumY;

		// auto leftRet = std::async(std::launch::async, &SchedulerService::consitencyCheckThread, this, out, leftSumY
		// 		, minTidx, minTidx+workLeft, ydivW, ymodW);
		auto fct = std::bind(&SchedulerService::consitencyCheckThread, this, out, system, std::ref(leftSumY)
				, minTidx, minTidx+workLeft, ydivW, ymodW);
		auto leftRet = std::async(std::launch::async, fct);

		ret = consitencyCheckThread(out, system, sumY, minTidx+workLeft, supTidx, ydivW, ymodW);
		ret = (ret & leftRet.get()); // wait for left thread

		for(int a = 0; a < sumY.size(); ++a)
			sumY[a] += leftSumY[a];
	}
	else //(minTidx == supTidx-1)
	{
		const int yminW = ydivW*(minTidx) + (ymodW > minTidx ? minTidx : ymodW );
		const int ysupW = yminW + ydivW + (ymodW > minTidx ? 1 : 0);

		sumY.resize(m_size.dimX(), 0);

		size_t cell = m_size.indexW(0, yminW);
		const size_t dimXW = m_size.dimXW();
		LatticePoint loc(0,0);
#define WITH_LOOP_CHECK
#ifdef WITH_LOOP_CHECK
		const SystemSizeCPU sizeCache(m_size);
#endif
		for(int yw = yminW; yw < ysupW; ++yw)
		{
			int sumX[BASIC_CELL_DIM_Y] {0,0,0,0};
			for(int xw = 0; xw < dimXW; ++xw, ++cell)
			{
				const unsigned int locSys = system[cell];
				for(int y = 0; y < BASIC_CELL_DIM_Y; ++y)
					for(int x = 0; x < BASIC_CELL_DIM_X; ++x)
					{
						loc.shift = SystemSize::shift(x,y);
						int loop = (loc(locSys) & 1) ? 1 : -1;
						sumX[y] += loop;
						const int gx = (xw<<L_BASIC_CELL_DIM_X) + x;
#ifdef WITH_LOOP_CHECK
						const int gy = (yw<<L_BASIC_CELL_DIM_Y) + y;
#endif
						{
							const int dy = (loc(locSys) & 2) ? 1 : -1;
							sumY[gx] += dy;
#ifdef WITH_LOOP_CHECK
							loop -= dy;
#endif
						}
#ifdef WITH_LOOP_CHECK
						loop += m_size.slopeY((gx+sizeCache.mDimX())&sizeCache.mDimX(), gy)(system) ? 1 : -1;
						loop -= m_size.slopeX(gx, (gy+sizeCache.mDimY())&sizeCache.mDimY())(system) ? 1 : -1;
						if(loop != 0)
						{
							ret = false;
							std::unique_lock<std::mutex> lock(out->second);
							out->first << "[MM][consitencyCheck] Found inconsistency in loop x= "
								<< gx << " , y= " << gy << " , offset= " << loop << '\n';
						}
#undef WITH_LOOP_CHECK
#endif
					}
			}
			for(int y = 0; y < BASIC_CELL_DIM_Y; ++y)
			{
				if(sumX[y] != 0)
				{
					ret = false;
					std::unique_lock<std::mutex> lock(out->second);
					out->first << "[MM][consitencyCheck] Found inconsistency in line y= " << (yw<<L_BASIC_CELL_DIM_Y)+y
						<< " , offset= " << sumX[y] << '\n';
				}
			}
		}
	}
	return ret;
}

bool Kpz::SchedulerService::consitencyCheck(bool full, std::ostream& out) const
{
	std::pair<std::ostream&, std::mutex> mout(std::piecewise_construct, std::forward_as_tuple(out), std::forward_as_tuple<>());
	std::vector<int> sumY;
	int nt = Correlator::nThreads();
	const int ydivW = m_size.dimYW()/nt;
	const int ymodW = m_size.dimYW()%nt;
	if(nt < 1)
		nt = 1;

	bool ret = consitencyCheckThread(&mout, m_system, sumY, 0, nt, ydivW, ymodW);

	for(int x = 0; x < sumY.size(); ++x)
	{
		if(sumY[x] != 0)
		{
			ret = false;
			out << "[MM][consitencyCheck] Found inconsistency in column x= " << x
				<< " , offset= " << sumY[x] << '\n';
		}
	}
	if(!ret)
		out << "[MM][consitencyCheck] Found inconsistency(s) in m_system\n";

	if(full)
	{
		for(auto corr : m_correlateTags)
		{
			std::fill(sumY.begin(), sumY.end(), 0);
			bool subret = consitencyCheckThread(&mout, corr->snapshot(), sumY, 0, nt,  ydivW, ymodW);
			for(int x = 0; x < sumY.size(); ++x)
			{
				if(sumY[x] != 0)
				{
					subret = false;
					out << "[MM][consitencyCheck] Found inconsistency in column x= " << x
						<< " , offset= " << sumY[x] << '\n';
				}
			}
			if(!subret)
				out << "[MM][consitencyCheck] Found inconsistency(s) in correlate tag mcs= " << corr->mcs() << '\n';
			ret = ret & subret;
		}
	}

	return ret;
}

std::ostream& Kpz::Roughness::printCompatible(std::ostream& o)
{
	return o << w2 << '\t' << h[0] << '\t' << h[1];
}

std::ostream& Kpz::Roughness::printMoments(std::ostream& o, unsigned int min, unsigned int sup)
{
	if(sup > H_ORDERS || min >= sup)
		throw std::range_error("[Roughness::printMoments] sup or min out of range.\n");
	o << h[min];
	for(int i = min+1; i < sup; ++i)
		o << '\t' << h[i];
	return o;
}
