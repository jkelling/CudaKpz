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

#include "schedulerCL.h"

#include <kmcRandom.h>

#include <unistd.h>

void Kpz::SchedulerCL::cleanupSize()
{
	clReleaseMemObject(d_System);
	clReleaseMemObject(d_Random);
	clReleaseProgram(m_program);
	clReleaseKernel(m_kpz);

	delete[] m_random;
}

Kpz::SchedulerCL::~SchedulerCL()
{
	//TODO: find out which objects have to be deleted here and how
	free(m_devices);
	clReleaseContext(m_context);

	cleanupSize();
}

void Kpz::SchedulerCL::pushSystem()
{
	cl_int err;
	err = clEnqueueWriteBuffer(m_commandQueue, d_System, CL_TRUE, 0, m_size.sizeW()*sizeof(unsigned int), m_system, 0, 0, 0);
	checkErr(err, "clEnqueueWriteBuffer(writeQueue)");
	m_deviceMcs = m_mcs;
}

void Kpz::SchedulerCL::popSystem()
{
	cl_int err;
	err = clEnqueueReadBuffer(m_commandQueue, d_System, CL_TRUE, 0, m_size.sizeW()*sizeof(unsigned int), m_system, 0, 0, 0);
	checkErr(err, "clEnqueueWriteBuffer(writeQueue)");
	m_mcs = m_deviceMcs;
}

void Kpz::SchedulerCL::mcs(int n)
{
	cl_int err;
	cl_event kernelStat;
	unsigned long long int N = n << L_MCS_DIV;
	for(; N > 0; --N)
	{
		const int xw = (int)(dsfmt_genrand_close_open(&m_dsfmt)*m_layout.blockDimXW());
		const int yw = (int)(dsfmt_genrand_close_open(&m_dsfmt)*m_layout.blockDimYW());
		err = clSetKernelArg(m_kpz, 2, sizeof(int), &xw);
		checkErr(err, "kernel.setArg(m_kpz, xw)");
		err = clSetKernelArg(m_kpz, 3, sizeof(int), &yw);
		checkErr(err, "kernel.setArg(m_kpz, yw)");
		err = clEnqueueNDRangeKernel( m_commandQueue, m_kpz, ThreadLayout::WORK_DIM, 0, m_layout.globalWorkSize(), m_layout.localWorkSize(), 0, NULL, &kernelStat);
		checkErr(err, "clEnqueueNDRangeKernel(m_kmc)");
		err = clWaitForEvents(1, &kernelStat);
		checkErr(err, "clWaitForEvents()");
		clReleaseEvent(kernelStat);
		++m_device;
	}
}

void Kpz::SchedulerCL::randomize()
{
#ifdef KMC_STRIDED_LCG_64
	m_random[0] = (unsigned long long)(dsfmt_genrand_close_open(&m_dsfmt)*((unsigned)-1))
		^ (((unsigned long long)(dsfmt_genrand_close_open(&m_dsfmt)*((unsigned)-1)))<<32);
	for (int a = 1; a<m_layout.getRandomNumberCount(); ++a)
	{
		m_random[a] = SLCGen(m_random[a-1]);
	}
#else
	for (int a = 0; a<m_layout.getRandomNumberCount(); ++a)
	{
		m_random[a] = (int)(dsfmt_genrand_close_open(&m_dsfmt)*(KMC_LCG_RAND_SUP))+1;
	}
#endif
	
	cl_int err;
	err = clEnqueueWriteBuffer(m_commandQueue, d_Random, CL_TRUE, 0, sizeof(KmcRandom_t)*m_layout.getRandomNumberCount(), m_random, 0, 0, 0);
	checkErr(err, "clEnqueueWriteBuffer(writeQueue)");
}

void Kpz::SchedulerCL::initOCL(int mode)
{
	m_mode = mode;
	cl_int err;

	cl_uint numPlatforms;
	
	err = clGetPlatformIDs(0, NULL, &numPlatforms);
	checkErr(err, "clGetPlatformIDs(numPlatforms)");

	cl_uint numDevices = 0;

	if (0 < numPlatforms) 
	{
		cl_platform_id* platforms = new cl_platform_id[numPlatforms];
		err = clGetPlatformIDs(numPlatforms, platforms, NULL);
		checkErr(err, "clGetPlatformIDs(platforms)");

		for (unsigned i = 0; i < numPlatforms; ++i) 
		{
			m_platform = platforms[i];

			switch(mode)
			{
				case MODE_GPU:
					err = clGetDeviceIDs(m_platform, CL_DEVICE_TYPE_GPU, 0, 0, (cl_uint*)&numDevices);
					break;
				case MODE_CPU:
				case MODE_CPU_DBG:
					err = clGetDeviceIDs(m_platform, CL_DEVICE_TYPE_CPU, 0, 0, (cl_uint*)&numDevices);
					break;
				default:
					std::cerr << "Unsupported mode.\n";
					exit(1);
			}
			if(numDevices > 0) // pick first device we find
				break;
		}
		delete[] platforms;
		if(numDevices == 0)
		{
			std::cerr << "No suitable device found.\n";
			exit(1);
		}
	}
	/*
	* If we could find our platform, use it. Otherwise pass a NULL and get whatever the
	* implementation thinks we should be using.
	*/
	{
		char name[100];
		err = clGetPlatformInfo(m_platform, CL_PLATFORM_VENDOR, 100, name, 0);
		checkErr(err, "clGetPlatformInfo(vendor)");
		cout() << "Vendor of choosen platform: " << name
			<< "\n-- number of GPUs: " << numDevices << '\n';
	}

	cl_context_properties cps[3] = {CL_CONTEXT_PLATFORM, (cl_context_properties)m_platform, 0};
	/* Use NULL for backward compatibility */
	cl_context_properties* cprops = (m_platform == 0) ? 0: cps;
	// End of platform layer
	
	m_devices = (cl_device_id*)malloc(numDevices * sizeof(cl_device_id));
	switch(mode)
	{
		case MODE_GPU:
			err = clGetDeviceIDs(m_platform, CL_DEVICE_TYPE_GPU, numDevices, m_devices, 0);
			break;
		case MODE_CPU:
		case MODE_CPU_DBG:
			err = clGetDeviceIDs(m_platform, CL_DEVICE_TYPE_CPU, numDevices, m_devices, 0);
			break;
	}
	checkErr(err, "clGetDeviceIDs");
	for(m_device = 0; m_device < numDevices; ++m_device)
	{
		m_context = clCreateContext(cprops, 1 /*numDevices*/, &m_devices[m_device], NULL, NULL, &err);
		//err = clGetDeviceInfo(m_devices[m_device], CL_DEVICE_AVAILABLE, sizeof(cl_bool), &av, 0);
		if(m_context)
			break;
	}
	if(m_device >= numDevices)
	{
		std::cerr << "FATAL ERROR: no GPU available.\n";
		checkErr(err, "clCreateContext()");
		exit(1);
	}

	switch(m_mode)
	{
		case MODE_CPU:
		case MODE_GPU:
			err = clGetDeviceInfo(m_devices[m_device], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(int), &m_blocks, 0);
			cout() << "-- using device " << m_device << " with " << m_blocks << " compute units\n";
			break;
		case MODE_CPU_DBG:
			m_blocks = 1;
			break;
	}
	size_t threads;
	err = clGetDeviceInfo(m_devices[m_device], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &threads, 0);
	m_code.setSkip(m_blocks, threads);

	m_code.loadPatchCode(std::string(getenv("KMCSVNDIR")) + "/kpz/kpz.cu");
	{
		std::ofstream file("raw.cl");
		file << m_code.raw();
		file.close();
	}
}

bool Kpz::SchedulerCL::Parser::paste(std::ostream& o, const std::string& directive) 
{
	if(directive.find("KMC_SLCG_SKIP") != std::string::npos)
	{
		o << "#define KMC_SLCG_A_SKIP (" << SLCGskipA(m_skip) << "lu)"
			<< "\n#define KMC_SLCG_C_SKIP (" << SLCGskipC(m_skip) << "lu)"
			<< "\n\n";
	}
	else
		return false;
	return true;
}

bool Kpz::SchedulerCL::setupOCL()
{
	size_t threads;
	cl_int err = clGetDeviceInfo(m_devices[m_device], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &threads, 0);
	size_t lmem;
	err = clGetDeviceInfo(m_devices[m_device], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(size_t), &lmem, 0);
	m_layout.setOptimal(threads, lmem, m_blocks, m_size);

	m_program = m_code.makeProgram(m_context, &m_devices[m_device], m_layout, m_size, m_kernelFlags);
	m_kpz = clCreateKernel(m_program, "kpzKernel", &err);
	checkErr(err, "clCreateKernel(kmcBlockDyn)");

	d_System = clCreateBuffer( m_context, CL_MEM_READ_WRITE, m_size.sizeW()*sizeof(unsigned int), 0, &err);
	checkErr(err, "d_System = clCreateBuffer()");
	m_random = new KmcRandom_t[m_layout.getRandomNumberCount()];
	d_Random = clCreateBuffer( m_context, CL_MEM_READ_WRITE, sizeof(KmcRandom_t)*m_layout.getRandomNumberCount(), 0, &err);
	checkErr(err, "d_Random = clCreateBuffer()");

	err = clSetKernelArg(m_kpz, 0, sizeof(d_System), &d_System);
	checkErr(err, "kernel.setArg(m_kmc, d_System)");
	err = clSetKernelArg(m_kpz, 3, sizeof(d_Random), &d_Random);
	checkErr(err, "kernel.setArg(m_kmc, d_Random)");

	m_commandQueue = clCreateCommandQueue(m_context, m_devices[m_device], 0, &err);
	checkErr(err, "clCreateCommandQueue()");
	return true;
}

bool Kpz::SchedulerCL::changedSizeEvent()
{
	if(!SchedulerService::changedSizeEvent())
		return false;
	cleanupSize();
	return setupOCL();
}

cl_program Kpz::SchedulerCL::Parser::makeProgram(cl_context context, cl_device_id* device, const ThreadLayout& layout, const SystemSize& size, const int flags)
{
	m_code.str("");

	//flags
	if(flags & KF_WAIT_IF_DEAD)
		m_code << "#define WAIT_IF_DEAD\n";
	if(flags & KF_OPTIMIZED_WAIT)
		m_code << "#define OPTIMIZED_WAIT\n";
	if(flags & KF_INNER_DC)
		m_code << "#define KPZ_INNER_DC\n";

	layout.printConst(m_code);
	size.printConst(m_code);
	m_code << m_raw;
	{
		std::ofstream file("code.cl");
		file << m_code.str();
		file.close();
	}
	return KmcCl::CuClParser::makeProgram(context, device);
}
