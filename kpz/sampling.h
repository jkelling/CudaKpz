/***************************************************************************
*   Copyright 2011 - 2016 Jeffrey Kelling <j.kelling@hzdr.de>
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

#include <memory>
#include <cmath>

namespace Kpz
{
	class FrontendArgs;
	class SchedulerService;

	class SamplingEvent
	{
		protected:
		unsigned int m_startTime;
		bool m_fCountAsMeasurement : 1;

		public:
			SamplingEvent(unsigned int startTime, bool fCountAsMeasurement = false)
				: m_startTime(startTime), m_fCountAsMeasurement(fCountAsMeasurement) {}

		virtual ~SamplingEvent() {}

		unsigned int startTime() const {return m_startTime;}
		bool fCountAsMeasurement() const {return m_fCountAsMeasurement;}

		virtual void handle(FrontendArgs& args, SchedulerService* scheduler) = 0;
	};

	class SamplingBarrier : public SamplingEvent
	{
		public:
			using SamplingEvent::SamplingEvent;

		virtual void handle(FrontendArgs& args, SchedulerService* scheduler);
	};

	class EndEvent : public SamplingEvent
	{
		public:
			using SamplingEvent::SamplingEvent;

		virtual void handle(FrontendArgs& args, SchedulerService* scheduler);
	};

	struct SamplingParam : public SamplingEvent
	{
		double measurementDensity, randomSamplingMean;
		unsigned int measurementInterval;


			SamplingParam(unsigned int startTime, double measurementDensity = .1, unsigned int measurementInterval = 0
					, double randomSamplingMean = NAN);

		bool isExp() const {return !std::isnan(measurementDensity);}
		bool isRandom() const {return !std::isnan(randomSamplingMean);}
		bool isLinear() const {return measurementInterval > 0;}

		bool setExp(double mD = .1);
		bool setRandom(double mean);
		bool setLinear(int i);

		template<class Scheduler>
		unsigned int nextEval(unsigned int mcs, Scheduler* scheduler) const {
			unsigned int ret;
			if(isLinear())
			{
				ret = mcs + measurementInterval;
			}
			else if(isExp())
			{
				ret = (unsigned int)ceil(exp(log(mcs+10)+measurementDensity));
			}
			else
			{
				ret = mcs + scheduler->genrand_close_open()*randomSamplingMean+randomSamplingMean;
			}
			return ret;
		}

		template<class Scheduler>
		unsigned int nextEvalFastForward(unsigned int mcs, Scheduler* scheduler) const {
			int ret = nextEval(mcs, scheduler);
			while(ret <= scheduler->mcs())
				ret = nextEval(ret, scheduler);
			return ret;
		}

		static std::shared_ptr<SamplingParam> makeExp(unsigned int startTime, double measurementDensity);
		static std::shared_ptr<SamplingParam> makeLinear(unsigned int startTime, unsigned int measurementInterval);
		static std::shared_ptr<SamplingParam> makeRandom(unsigned int startTime, double mean);

		virtual void handle(FrontendArgs& args, SchedulerService* scheduler);
	};

	class CorrelateTag : public SamplingEvent
	{
		unsigned int m_stopTime;

		public:

			CorrelateTag(unsigned int startTime, unsigned int stopTime = 0, bool fCountAsMeasurement = false)
				: SamplingEvent(startTime, fCountAsMeasurement), m_stopTime(stopTime)
			{}

		unsigned int stopTime() const {return m_stopTime;}

		virtual void handle(FrontendArgs& args, SchedulerService* scheduler);
	};
}
