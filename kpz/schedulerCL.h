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

#ifndef KPZ_SCHEDULER_CL_H
#define KPZ_SCHEDULER_CL_H

#include "schedulerService.h"
#include "systemSize.h"
#include "threadLayout.h"

#include <opencl/util.h>

#include <string>

namespace Kpz
{
class SchedulerCL : public SchedulerService
{
	class Parser : public KmcCl::CuClParser
	{
		std::string m_raw;
		int m_skip;
		protected:
		bool paste(std::ostream& o, const std::string& directive);
		public:
		int loadPatchCode(const std::string &file) {
			int n = KmcCl::CuClParser::loadPatchCode(file);
			m_raw = m_code.str();
			return n;
		}
		cl_program makeProgram(cl_context context, cl_device_id* device, const ThreadLayout& layout, const SystemSize& size, const int flags);

		const std::string& raw() const {return m_raw;}
		void setSkip(int maxBlocks, int maxThreads = 1024) {m_skip = maxBlocks*maxThreads;}
	};

	ThreadLayout m_layout;
	int m_blocks;

	cl_platform_id m_platform;
	cl_device_id* m_devices;
	int m_device;
	cl_context m_context;
	cl_program m_program;
	cl_kernel m_kpz;
	cl_mem d_System;
	cl_mem d_Random;
	cl_command_queue m_commandQueue;
	KmcRandom_t* m_random;

	Parser m_code;
	int m_mode, m_kernelFlags;

	void initOCL(int mode); // init opencl handles and device mem
	bool setupOCL();
	void cleanupSize();

	// This function is only called when the systemsize is changed after the class
	// has been constructed. As long as the load functions are not implemented this
	// is unlikly to ever occour.
	virtual bool changedSizeEvent(); // call base class::setSize(lx,ly) first

	public:
	static const int MODE_GPU = 1;
	static const int MODE_CPU_DBG = 2;
	static const int MODE_CPU = 3;
	static const int MODE_DEFAULT = MODE_GPU;

	static const int KF_WAIT_IF_DEAD = 1;
	static const int KF_OPTIMIZED_WAIT = 2;
	static const int KF_INNER_DC = 4;
	static const int KF_DEFAULT = 0
#ifdef WAIT_IF_DEAD
		| KF_WAIT_IF_DEAD
#endif
#ifdef OPTIMIZED_WAIT
		| KF_OPTIMIZED_WAIT
#endif
#ifdef KPZ_INNER_DC
		| KF_INNER_DC
#endif
		;

		SchedulerCL(int lx, int ly, int kernelFlags = KF_DEFAULT, int mode = MODE_DEFAULT)
			: SchedulerService(lx, ly), m_kernelFlags(kernelFlags)
		{initOCL(mode);setupOCL();}
		SchedulerCL(const SystemSize& size, int kernelFlags = KF_DEFAULT, int mode = MODE_DEFAULT)
			: SchedulerService(size), m_kernelFlags(kernelFlags)
		{initOCL(mode);setupOCL();}
		~SchedulerCL(); // cleaup ocl

	void pushSystem(); // put system on device
	void popSystem(); // get system form device
	void mcs(int n); // do n mcs
	void randomize(); // init on device rng
};

}

#endif
