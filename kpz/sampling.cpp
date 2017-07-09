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

#include "sampling.h"

#include "kpzFrontend.h"
#include "schedulerService.h"

Kpz::SamplingParam::SamplingParam(unsigned int startTime, double measurementDensity, unsigned int measurementInterval
		, double randomSamplingMean)
	: SamplingEvent(startTime, true)
		, measurementDensity(measurementDensity), randomSamplingMean(randomSamplingMean)
		, measurementInterval(measurementInterval)
{}

bool Kpz::SamplingParam::setExp(double mD)
{
	if(std::isnan(mD))
		return false;
	measurementDensity = mD;
	measurementInterval = 0;
	randomSamplingMean = NAN;
	return true;
}

bool Kpz::SamplingParam::setRandom(double mean)
{
	if(std::isnan(mean))
		return false;
	measurementDensity = NAN;
	measurementInterval = 0;
	randomSamplingMean = mean;
	return true;
}

bool Kpz::SamplingParam::setLinear(int i)
{
	if(i == 0)
		return false;
	measurementDensity = NAN;
	measurementInterval = i;
	randomSamplingMean = NAN;
	return true;
}

std::shared_ptr<Kpz::SamplingParam> Kpz::SamplingParam::makeExp(unsigned int startTime, double measurementDensity)
{
	auto ret = std::make_shared<SamplingParam>(startTime);
	ret->setExp(measurementDensity);
	return ret;
}

std::shared_ptr<Kpz::SamplingParam> Kpz::SamplingParam::makeLinear(unsigned int startTime, unsigned int measurementInterval)
{
	auto ret = std::make_shared<SamplingParam>(startTime);
	ret->setLinear(measurementInterval);
	return ret;
}

std::shared_ptr<Kpz::SamplingParam> Kpz::SamplingParam::makeRandom(unsigned int startTime, double mean)
{
	auto ret = std::make_shared<SamplingParam>(startTime);
	ret->setRandom(mean);
	return ret;
}

void Kpz::SamplingBarrier::handle(FrontendArgs& args, SchedulerService* scheduler)
{
	args.setSamplingSuspended(false);
}

void Kpz::EndEvent::handle(FrontendArgs& args, SchedulerService* scheduler)
{
	args.setFinished(true);
}

void Kpz::SamplingParam::handle(FrontendArgs& args, SchedulerService* scheduler)
{
	args.setCurrentSampling(std::make_shared<SamplingParam>(*this));
}

void Kpz::CorrelateTag::handle(FrontendArgs& args, SchedulerService* scheduler)
{
	scheduler->correlateTag(m_stopTime);
}
