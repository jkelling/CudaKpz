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

#ifndef KPZ_CORRELATOR_H
#define KPZ_CORRELATOR_H

#include "schedulerService.h"

#include <iostream>

namespace splash
{
	class DataCollector;
}

namespace Kpz
{
	class Correlator
	{
		protected:
		struct CorrelateData;
		static void* correlateThreadCO(void* cd);
		static void* roughnessThreadCO(void* cd);
		static void* correlateThreadCBBit(void* cd);
		static void* roughnessThreadCBBit(void* cd);

		SchedulerService* m_scheduler;
		unsigned int* m_snapshot;
		double m_h;
		double m_ch, m_cs;
		unsigned int m_mcs, m_stop;

			Correlator(const Correlator&); //no copy

		static unsigned int s_nThreads;

		public:

			Correlator(SchedulerService* scheduler, unsigned int stop = std::numeric_limits<unsigned int>::max());
			virtual ~Correlator() {delete[] m_snapshot;}

		virtual Roughness set(unsigned int mcs);
		virtual Roughness correlate();
		virtual Roughness correlateThreaded();
#ifdef MODE_CUDA
		virtual Roughness correlateCuda(unsigned int* d_System, int nBlocks);
#else
		virtual Roughness correlateCuda(unsigned int* d_System, int nBlocks) {exit(1);}
#endif
		virtual Roughness roughness();
		virtual Roughness roughnessThreaded();
		void save(splash::DataCollector* data, const std::string& prefix);
	
		double h() const {return m_h;}
		double ch() const {return m_ch;}
		double cs() const {return m_cs;}
		unsigned int mcs() const {return m_mcs;}
		unsigned int stop() const {return m_stop;}
		const unsigned int* snapshot() const {return m_snapshot;}
		
		static Correlator* fromH5(splash::DataCollector* data, int id, SchedulerService* scheduler
				, const std::string& prefix);

		static unsigned int nThreads() {return s_nThreads;}
		static void setNThreads(unsigned int n);
	};
}

std::ostream& operator<< (std::ostream& o, const Kpz::Correlator& c);

#endif
