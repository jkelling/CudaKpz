/***************************************************************************
*   Copyright 2014 - 2014 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef KPZ_GPU_TINY_MT_H
#define KPZ_GPU_TINY_MT_H

#define DSFMT_MEXP 19937
#include <dSFMT/dSFMT.h>

#include <string>

namespace splash {
	class DataCollector;
}

class TinyMT32
{
	unsigned int* d_Random, *d_mat1, *d_mat2, *d_tmat;
	size_t m_generatorCount;
	dsfmt_t* m_dsfmt;
	int m_blocks, m_threads;

	void initCUDA();

	public:
		TinyMT32(size_t generatorCount, dsfmt_t* dsfmt)
			: m_generatorCount(generatorCount), m_dsfmt(dsfmt) {initCUDA();}
		~TinyMT32();

	class Device;

	void randomize();

	void* getDRandom(int numbersPerGen) {return d_Random;}
	size_t generatorCount() const {return m_generatorCount;}

	void writeH5(splash::DataCollector* data, int id, const std::string& prefix = "");
	bool readH5(splash::DataCollector* data, int id, const std::string& prefix = "");
};
#endif
