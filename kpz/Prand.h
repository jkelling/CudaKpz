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

#ifndef KPZ_PRAND_H
#define KPZ_PRAND_H

#define DSFMT_MEXP 19937
#include <dSFMT/dSFMT.h>

const size_t PRNG_L_MAX_BUFFER_SIZE = 31;

#ifndef __CUDACC__
class mt19937_state;
#else
#include <mt19937.h>
#endif

#include <iostream>

/*! [Comput. Phys. Commun. 185(2014)1343]( http://cpc.cs.qub.ac.uk/summaries/AESB_v1_0.html )
 */
class PrandMT
{
	unsigned int* d_Random;
	size_t m_generatorCount, m_bufferSize, m_used;
	mt19937_state* m_state;
	dsfmt_t* m_dsfmt;

	void initCUDA();
	bool fillBuffer(int numbers);

	public:
		PrandMT(size_t generatorCount, dsfmt_t* dsfmt);
		~PrandMT();

	class Device;

	void randomize();

	inline void* getDRandom(int numbersPerGen) {
		const int numbers = numbersPerGen*m_generatorCount;
		if(m_used+numbers > m_bufferSize)
			if(!fillBuffer(numbers))
				return 0;
		void* ret = &d_Random[m_used];
		m_used += numbers;
		return ret;
	}
};

#endif
