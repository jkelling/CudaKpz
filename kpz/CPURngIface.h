/***************************************************************************
*   Copyright 2014 - 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef KPZ_PARALLEL_CPU_RNG_IFACE_H
#define KPZ_PARALLEL_CPU_RNG_IFACE_H

#include <dSFMT/dSFMT.h>
#include <tinyMT/tinymt32.h>

#include <iostream>
#include <climits>

namespace splash {
	class DataCollector;
}

namespace Kpz
{
namespace RngIface
{
	/*! \brief Provides a generic interfaces to dSFMT.
	 *
	 * \warning Do not use with multiple threads, parallel sequences may be
	 * correlated.
	 * \todo adjust to new interface required by SchedulerCPUBaseThreadedRng
	 */
	class DSFMT
	{
		dsfmt_t m_dsfmt;

		public:
			
		// general interface
		void seed(unsigned long long seed) {
			dsfmt_init_gen_rand(&m_dsfmt, seed);
			for( int a = 0; a < 10000; ++a) // arbitray large number 
			{
				dsfmt_genrand_close_open(&m_dsfmt);
			}
		}

		double close_open() {return dsfmt_genrand_close_open(&m_dsfmt);}

		static DSFMT* create(int n, dsfmt_t* dsfmt);

		// special
		dsfmt_t* getDSFMT() {return &m_dsfmt;}
	};

	/*! \brief Provides a generic interfaces to tinyMT.
	 *
	 */
	class TinyMT32
	{
		TINYMT32_T m_tmt;

		public:

			TinyMT32(int index) {setMat(index);}
			TinyMT32() {}

		void setMat(int index);

		// general interface
		void seed(unsigned int seed) {
			tinymt32_init(&m_tmt, seed);
		}

		unsigned int random() {return tinymt32_generate_uint32(&m_tmt);}
		float randomFloat() {return tinymt32_generate_float(&m_tmt);}

		template<class Container>
		static void initialize(Container& arr) {
			int i = 0;
			for(auto a = arr.begin(); a != arr.end(); ++a, ++i)
				a->setMat(i);
			std::cout << "TinyMT32: Using parametersets: 0 .. " << i-1 << '\n';
		}

		template<class Container>
		static void randomize(Container& arr, dsfmt_t& dsfmt) {
			for(auto a = arr.begin(); a != arr.end(); ++a)
				a->seed(dsfmt_genrand_close_open(&dsfmt)*UINT_MAX);
		}

		template<class Container>
		static void writeH5(splash::DataCollector* data, int id, Container& arr, const std::string& prefix = "") {}
		template<class Container>
		static bool readH5(splash::DataCollector* data, int id, Container& arr, const std::string& prefix = "")
		{return false;}
	};
}
}

#endif
