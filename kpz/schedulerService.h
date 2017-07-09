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

#ifndef KPZ_SCHEDULER_SERVICE_H
#define KPZ_SCHEDULER_SERVICE_H

#define DSFMT_MEXP 19937
#include <dSFMT/dSFMT.h>

#include <iostream>
#include <pthread.h>
#include <list>
#include <limits>
#include <future>
#include <mutex>
#include <vector>

#include "systemSize.h"

namespace splash
{
	class DataCollector;
}

namespace MPI
{
	class Comm;
}

namespace Kpz
{
	class Correlator;

	struct Roughness {
		static const int H_ORDERS = 7;
		long double w2, h[H_ORDERS];

		const Roughness& operator+=(const Roughness& other);
		const Roughness& operator/=(double div);
		inline void addH(double ph) {
			h[0] += ph;
			double acc = ph;
			for(int a = 1; a < H_ORDERS; ++a)
				h[a] += (acc*=ph);
		}
		void normAndSetW2(unsigned long long N);

		std::ostream& printCompatible(std::ostream& o);
		std::ostream& printMoments(std::ostream& o, unsigned int min = 0, unsigned int sup = H_ORDERS);

			Roughness() : w2(0.) {for(auto& a : h) a=0.;}
			Roughness(const Roughness* r, int n);
	};

	class SchedulerService
	{
		public:
		enum Encoding {ENC_LOCAL, ENC_CBBIT};
		protected:

		struct CallMcsArgs {
			SchedulerService* scheduler;
			void (SchedulerService::*fkt) (int);
			int n;
		} m_mcsArgs;
		static void* callMcs(void* args);
		void* (*m_callMcs)(void*);
		inline void mcsAsync(int n, void (SchedulerService::*fkt)(int))
		{
			joinMcs();
			m_mcsArgs.scheduler = this;
			m_mcsArgs.fkt = fkt;
			m_mcsArgs.n = n;
			pthread_create(&m_mcsThread, NULL, m_callMcs, (void*) &m_mcsArgs);
		}

		Encoding m_encoding;
		unsigned int* m_system, *m_disorder;
		dsfmt_t m_dsfmt;
		SystemSize m_size;

		unsigned int m_deviceMcs;
		pthread_t m_roughnessThread, m_mcsThread;
		Roughness m_roughness;
		static void* callRoughness(void* scheduler);
		static void* threadCorrelateTag(void* scheduler);
		static void* threadCorrelate(void* scheduler);
		static void* threadCorrelateLinesTag(void* scheduler);
		bool m_silentRoughness;
		int m_outPrecision;
		bool consitencyCheckThread(std::pair<std::ostream&, std::mutex>* out, const unsigned int* system, std::vector<int>& sumY
				, int minTidx, int supTidx, int ydivW, int ymodW) const;

		bool m_fAnalysisOnGPU : 1;

		std::ostream* out;

		// simulation params
		double m_disorderP[4], m_disorderQ[4];
		static const int RNG_SYNC_PACAGE = 3;
		struct RngPackage {
			double m[RNG_SYNC_PACAGE];

				RngPackage(dsfmt_t& dsfmt);
			const double& operator[] (int i) {return m[i];}
		};

		unsigned int m_mcs, m_stopCorr;
		std::list<Correlator*> m_correlateTags;

		void reSeed(int seed = 0);
		void initURandom();
		bool setSize();
		virtual bool changedSizeEvent() {return true;}

		bool loadXYZ(const char* file);
		bool loadBit(const char* file);

		public:

			SchedulerService(int lx, int ly, Encoding enc = ENC_LOCAL);
			SchedulerService(const SystemSize& size, Encoding enc = ENC_LOCAL);
			SchedulerService(SchedulerService &&other);
			virtual ~SchedulerService();

		virtual void initHomogeneous(Encoding encoding = ENC_LOCAL);
		void nullSystem();
		virtual Roughness roughness();
		Roughness roughnessLines();
		bool consitencyCheck(bool full = true, std::ostream& out = std::cout) const;

		void roughnessAsync();
		virtual void correlateTag(unsigned int stop = std::numeric_limits<unsigned int>::max());
		virtual void correlate();
		void setSilentRoughness(bool b = true) {m_silentRoughness = b;}
		void setOutPrecision(int outPrecision = 12) {m_outPrecision = outPrecision;}
		void setAnalysisOnGPU(bool b = true) {m_fAnalysisOnGPU = b;}
		const Roughness& getRoughness() {
			joinRoughness();
			return m_roughness;
		}
		const std::list<Correlator*>& correlateTags() {
			joinRoughness();
			return m_correlateTags;
		}

		inline void joinRoughness() {
			if(m_roughnessThread)
			{
				pthread_join(m_roughnessThread, NULL);
				m_roughnessThread = 0;
			}
		}
		inline void joinMcs() {
			if(m_mcsThread)
			{
				pthread_join(m_mcsThread, NULL);
				m_mcsThread = 0;
			}
		}
		inline void join() {
			joinRoughness();
			joinMcs();
		}

		bool setSize(int lx, int ly);
		bool initDisorder2();
		virtual bool setPQ(double p, double q, unsigned int n = 0) {
			if(n >= 4)
				return false;
			m_disorderP[n] = p;
			m_disorderQ[n] = q;
			return true;
		}
		double rateP(unsigned int n = 0) const {return m_disorderP[n];}
		double rateQ(unsigned int n = 0) const {return m_disorderQ[n];}
		bool generateDisorder2(double fill);

		virtual void pushSystem();
		virtual void popSystem();
		virtual void mcs(int n) {}
		void mcsAsync(int n) {mcsAsync(n, &SchedulerService::mcs);}
		/*! MCS with disorder
		 */
		virtual void mcsPQ(int n) {std::cerr << "[BUG] call to SchedulerService::mcsPQ()\n";}
		virtual void mcsPQAsync(int n) {mcsAsync(n, &SchedulerService::mcsPQ);}
		virtual void mcsDisorder2(int n) {std::cerr << "[BUG] call to SchedulerService::mcsDisorder2()\n";}
		virtual void mcsDisorder2Async(int n) {mcsAsync(n, &SchedulerService::mcsDisorder2);}

		virtual void mcsSyncRng(int n) {std::cerr << "[BUG] call to SchedulerService::mcsSyncRng()\n";}
		virtual void mcsSyncRngAsync(int n) {mcsAsync(n, &SchedulerService::mcsSyncRng);}
		virtual void mcsPQSyncRng(int n) {std::cerr << "[BUG] call to SchedulerService::mcsPQSyncRng()\n";}
		virtual void mcsPQSyncRngAsync(int n) {mcsAsync(n, &SchedulerService::mcsPQSyncRng);}
		virtual void mcsDisorder2SyncRng(int n) {std::cerr << "[BUG] call to SchedulerService::mcsDisorder2SyncRng()\n";}
		virtual void mcsDisorder2SyncRngAsync(int n) {mcsAsync(n, &SchedulerService::mcsDisorder2SyncRng);}

		virtual void mcsNop(int n) {m_deviceMcs += n;}

		virtual void randomize();
		
		bool saveBit(const char* file, const char* meta);
		bool saveXYZ(const char* file);
		bool load(const char* file);

		void supplySeed(int std, int dsfmt);
		std::ostream& cout() const {return *out;}
		virtual void setOut (std::ostream& o) {out = &o;}

		const SystemSize& size() const {return m_size;}
		unsigned int* system() {return m_system;}
		const unsigned int* system() const {return m_system;}
		unsigned int* disorder() {return m_disorder;}
		const unsigned int* disorder() const {return m_disorder;}

		const dsfmt_t& dsfmt() const {return m_dsfmt;}
		void copyDsfmt(const dsfmt_t& dsfmt) {m_dsfmt = dsfmt;}
		int syncRNGMPI(MPI::Comm& comm);
		double genrand_close_open() {return dsfmt_genrand_close_open(&m_dsfmt);}

		bool copySystem(const SchedulerService* s);
		bool copyDisorder(const SchedulerService* s);

		virtual double mcsScale() const {return 1.;}

		unsigned int mcs() const {return m_mcs;}
		virtual void setMcs(unsigned int mcs);
		unsigned int deviceMcs() const {return m_deviceMcs;}
		Encoding encoding() const {return m_encoding;}

		splash::DataCollector* writeH5(const char* file, const std::string& prefix = "");
		splash::DataCollector* readH5(const char* file, const std::string& prefix = "");
		static splash::DataCollector* readH5open(const char* file);
		virtual bool collectH5(splash::DataCollector* data, const std::string& prefix = "");
		virtual bool retrieveH5(splash::DataCollector* data, const std::string& prefix = "");
		static int mcsH5(splash::DataCollector* data);
		static void closeH5(splash::DataCollector* data);
	};
}
std::ostream& operator<<(std::ostream& o, const Kpz::Roughness& r);

#endif
