/***************************************************************************
*   Copyright 2014 - 2016 Jeffrey Kelling <j.kelling@hzdr.de>
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

#define DSFMT_MEXP 19937
#include "../dSFMT/dSFMT.h"

#include "../Random123/threefry.h"

#include <vector>
#include <string>

#ifdef USE_LIB_SPLASH
namespace splash {
	class DataCollector;
}
#endif
namespace H5 {
	class CommonFG;
}

/*! \warning generatorCount will be rounded up to be a multiple of the warp size.
 */
class Random123_Threefry
{
	public:
	// the default typedef in r123 has 20 rounds
	typedef	r123::Threefry2x32_R<20> Threefry;

	private:

	unsigned int* d_Random;
	size_t m_generatorCount;
	dsfmt_t* m_dsfmt;
	Threefry::key_type m_key;

	void initCUDA();

	std::string mkH5DsName(const char* prefix = 0);

	public:

		Random123_Threefry(size_t generatorCount, dsfmt_t* dsfmt = 0)
			: m_generatorCount(generatorCount), m_dsfmt(dsfmt) {initCUDA();}
		~Random123_Threefry();

	class Device;

	void randomize();
	// void randomize(const std::vector<unsigned int>& seeds);
	/*! Ensure that at least \p ngenmin generators are available.
	 * \return whether m_generatorCount sufficed without adjustment. */
	bool minGenerators(unsigned int ngenmin);

	void* getDRandom(int numbersPerGen) {return d_Random;}
	size_t generatorCount() const {return m_generatorCount;}

	/*! Call this function in each compilaton unit where a kernel uses the generator. */
	void initDeviceConst();

// #ifdef USE_LIB_SPLASH
// 	void writeH5(splash::DataCollector* data, int id, const char* prefix = 0);
// 	bool readH5(splash::DataCollector* data, int id, const char* prefix = 0);
// #endif
// 	void writeH5(H5::CommonFG& dest, const char* prefix = 0);
// 	void readH5(H5::CommonFG& src, const char* prefix = 0);
};
