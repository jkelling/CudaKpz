/***************************************************************************
*   Copyright 2015 - 2015 Jeffrey Kelling <j.kelling@hzdr.de>
*                  Helmholtz-Zentrum Dresden-Rossendorf
*                  Institute of Ion Beam Physics and Materials Research
*
*   This program is free software; you can redistribute it and/or
*   modify it under the terms of the GNU General Public
*   License as published by the Free Software Foundation; either
*   version 2 of the License, or (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program; if not, write to the Free Software
*   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
***************************************************************************/

#include <iomanip>
#include <climits>

template<class Rng>
void Kpz::SchedulerCPUCBBit<Rng>::Func_MCSLocalBase::pSyncRngUpdateVL5(int yMin, int ySup)
{
	// compiler error in gcc 4.8.2
#if 0
	static const int L_VECTOR_SIZE = 5; // AVX
	typedef unsigned long long uintVec __attribute__ ((vector_size (1<<L_VECTOR_SIZE)));
	uintVec *system = (uintVec*)m_this->system();
	// auto printUintVec = [](const uintVec& v, std::ostream& o = std::cerr) -> std::ostream&  {
	// 	for (int a = 0; a < 4; ++a)
	// 		o << v[a] << ' ';
	// 	return o;
	// };

	const size_t quadrant = m_this->size().sizeW()>>L_VECTOR_SIZE;
	// std::cerr << "loop xw< " << (m_this->m_size.dimX()>>(L_VECTOR_SIZE+3+1)) << " yMin= " << yMin << " ySup= " << ySup
	// 	<< " quadrant= " << quadrant << " systemptr= " << system << '\n';
	const int signum = m_this->m_signum;
	for(int y = yMin; y < ySup; ++y)
	{
		const bool parity = signum==(y&1);
		// xmax divided by another 2 for sublattice splitting
		for(int xw = 0; xw < (m_this->m_size.dimX()>>(L_VECTOR_SIZE+3+1)); ++xw)
		{
			const size_t offset = xw+(y<<(m_this->size().lDimX()-(L_VECTOR_SIZE+4)));
			// std::cerr << "xw,y: " << xw << ',' << y << " offset= " << offset << " signum= " << signum << '\n';
			size_t index = ((quadrant<<1)&(0-signum))+offset;
			// std::cerr << "thisData* index: " << index << ' ' << index + quadrant << '\n';
			uintVec thisDataX = system[index];
			uintVec thisDataY = system[index+quadrant];
			// printUintVec(thisDataX) <<  " =thisDataX\n";
			// printUintVec(thisDataY) <<  " =thisDataY\n";
			index = ((quadrant<<1)&(~(0-signum)));
			const size_t yIndexY = index + quadrant
				+ xw+ (((y+1)&(m_this->size().dimY()-1))<<(m_this->size().lDimX()-(L_VECTOR_SIZE+4)));
			// std::cerr << "*Data* index: " << index << ' ' << yIndexY << '\n';
			index += offset;
			uintVec yDataY = system[yIndexY];
			uintVec xDataX = system[index];
			// printUintVec(xDataX) <<  " =xDataX\n";
			// printUintVec(yDataY) <<  " =yDataY\n";
			// std::cerr << "*Data* index: " << index << ' ' << yIndexY << '\n';

			if(!parity)
			{
				// std::cerr << std::hex;
				index += (xw+1 == 1<<(m_this->size().lDimX()-(L_VECTOR_SIZE+3+1))) ? (-xw) : 1;
				auto tmp = xDataX << 63; // preserve LSBs of elements
				tmp[0] = ((*(unsigned long long*)(system + index))&1)<<63; // get LSB of next vector
				static const uintVec SMASK = {1,2,3, 0};
				tmp = __builtin_shuffle(tmp, SMASK); // rotate ints
				xDataX >>= 1;
				static const unsigned long long MSB = 1ull<<63;
				static const uintVec mask = {MSB,MSB,MSB,MSB};
				xDataX &= ~mask; // make sure MSBs are 0
				xDataX ^= tmp;
				// printUintVec(xDataX) << " =xDataX shifted + next LSB\n";
			}
			
			thisDataX = ~(thisDataX | thisDataY); // both x and y of this unset
			thisDataX &= xDataX & yDataY; // ... and both xx and yy set

			{
				//\todo use TinyMT64
				const uintVec rnd {
					((unsigned long long)m_rng->random()<<32)^m_rng->random(),
					((unsigned long long)m_rng->random()<<32)^m_rng->random(),
					((unsigned long long)m_rng->random()<<32)^m_rng->random(),
					((unsigned long long)m_rng->random()<<32)^m_rng->random(),
				};
				thisDataX &= rnd; // flip mask
			}

			xDataX = thisDataX;
			if(!parity)
			{
				// use index as computed above for !parity case
				auto tmp = xDataX >> 63; // preserve MSBs of elements
				static const uintVec SMASK = {3, 0,1,2};
				tmp = __builtin_shuffle(tmp, SMASK); // rotate ints back
				(*(unsigned long long*)(system + index)) ^= tmp[0]&1; // attempt to flip LSB of next vector
				tmp[0] = 0; // do not touch LSB of xDataX
				xDataX <<= 1;
				xDataX ^= tmp;
			}
			index = ((quadrant<<1)&(~(0-signum)))+offset;
			system[index] ^= xDataX;
			index = ((quadrant<<1)&(0-signum))+offset;
			system[index] ^= thisDataX;
			system[index+quadrant] = thisDataY ^ thisDataX;
			system[yIndexY] = yDataY ^ thisDataX;
		}
	}
#endif
}

template<class Rng>
void Kpz::SchedulerCPUCBBit<Rng>::Func_MCSLocalBase::pSyncRngUpdateVL4(int yMin, int ySup)
{
	// compiler error in gcc 4.8.2
#if 0
	static const int L_VECTOR_SIZE = 4; // SSE
	typedef unsigned int uintVec __attribute__ ((vector_size (1<<L_VECTOR_SIZE)));
	uintVec *system = (uintVec*)m_this->system();

	const size_t quadrant = m_this->size().sizeW()>>L_VECTOR_SIZE;
	const int signum = m_this->m_signum;
	for(int y = yMin; y < ySup; ++y)
	{
		const bool parity = signum==(y&1);
		// xmax divided by another 2 for sublattice splitting
		for(int xw = 0; xw < (m_this->m_size.dimX()>>(L_VECTOR_SIZE+3+1)); ++xw)
		{
			const size_t offset = xw+(y<<(m_this->size().lDimX()-(L_VECTOR_SIZE+4)));
			size_t index = ((quadrant<<1)&(0-signum))+offset;
			uintVec thisDataX = system[index];
			uintVec thisDataY = system[index+quadrant];
			index = ((quadrant<<1)&(~(0-signum)));
			const size_t yIndexY = index + quadrant
				+ xw+ (((y+1)&(m_this->size().dimY()-1))<<(m_this->size().lDimX()-(L_VECTOR_SIZE+4)));
			index += offset;
			uintVec yDataY = system[yIndexY];
			uintVec xDataX = system[index];

			if(!parity)
			{
				index += (xw+1 == 1<<(m_this->size().lDimX()-(L_VECTOR_SIZE+3+1))) ? (-xw) : 1;
				auto tmp = xDataX << 31; // preserve LSBs of elements
				tmp[0] = ((*(unsigned int*)(system + index))&1)<<31; // get LSB of next vector
				static const uintVec SMASK = {1,2,3, 0};
				tmp = __builtin_shuffle(tmp, SMASK); // rotate ints
				xDataX >>= 1;
				static const unsigned long long MSB = 1ull<<31;
				static const uintVec mask = {MSB,MSB,MSB,MSB};
				xDataX &= ~mask; // make sure MSBs are 0
				xDataX ^= tmp;
			}

			thisDataX = ~(thisDataX | thisDataY); // both x and y of this unset
			thisDataX &= xDataX & yDataY; // ... and both xx and yy set

			{
				//\todo use TinyMT64
				const uintVec rnd {
					m_rng->random(),
					m_rng->random(),
					m_rng->random(),
					m_rng->random(),
				};
				thisDataX &= rnd; // flip mask
			}

			xDataX = thisDataX;
			if(!parity)
			{
				// use index as computed above for !parity case
				auto tmp = xDataX >> 31; // preserve MSBs of elements
				static const uintVec SMASK = {3, 0,1,2};
				tmp = __builtin_shuffle(tmp, SMASK); // rotate ints back
				(*(unsigned int*)(system + index)) ^= tmp[0]&1; // attempt to flip LSB of next vector
				tmp[0] = 0; // do not touch LSB of xDataX
				xDataX <<= 1;
				xDataX ^= tmp;
			}
			index = ((quadrant<<1)&(~(0-signum)))+offset;
			system[index] ^= xDataX;
			index = ((quadrant<<1)&(0-signum))+offset;
			system[index] ^= thisDataX;
			system[index+quadrant] = thisDataY ^ thisDataX;
			system[yIndexY] = yDataY ^ thisDataX;
		}
	}
#endif
}

template<class Rng>
void Kpz::SchedulerCPUCBBit<Rng>::Func_MCSLocalBase::pSyncRngUpdateVL3(int yMin, int ySup)
{
	static const int L_VECTOR_SIZE = 3; // 64-bit
	typedef unsigned long long uintVec;
	uintVec *system = (uintVec*)m_this->stash();

	const size_t quadrant = m_this->size().sizeW()>>L_VECTOR_SIZE;
	const int signum = m_this->m_signum;
	for(int y = yMin; y < ySup; ++y)
	{
		const bool parity = signum==(y&1);
		// xmax divided by another 2 for sublattice splitting
		for(int xw = 0; xw < (m_this->m_size.dimX()>>(L_VECTOR_SIZE+3+1)); ++xw)
		{
			const size_t offset = xw+(y<<(m_this->size().lDimX()-(L_VECTOR_SIZE+4)));
			size_t index = ((quadrant<<1)&(0-signum))+offset;
			uintVec thisDataX = system[index];
			uintVec thisDataY = system[index+quadrant];
			index = ((quadrant<<1)&(~(0-signum)));
			const size_t yIndexY = index + quadrant
				+ xw+ (((y+1)&(m_this->size().dimY()-1))<<(m_this->size().lDimX()-(L_VECTOR_SIZE+4)));
			uintVec yDataY = system[yIndexY];
			index += offset;
			uintVec xDataX = system[index];

			if(!parity)
			{
				index += (xw+1 == 1<<(m_this->size().lDimX()-(L_VECTOR_SIZE+3+1))) ? (-xw) : 1;
				xDataX >>= 1;
				static const unsigned long long MSB = 1ull<<63;
				xDataX &= ~MSB; // delete MSB
				xDataX |= ((system[index]<<63)&MSB); // set MSB to LSB of next vector
			}
			
			const auto thisDataXBak = thisDataX;
			thisDataX = ~(thisDataX | thisDataY); // both x and y of this unset
			thisDataX &= xDataX & yDataY; // ... and both xx and yy set

			{
				//\todo use TinyMT64
				const uintVec rnd =	(((unsigned long long)m_rng->random()<<32)^m_rng->random());
				thisDataX &= rnd; // flip mask
			}

			xDataX = thisDataX;
			if(!parity)
			{
				// use index as computed above for !parity case
				system[index] ^= (thisDataX>>63)&1; // attempt flip of LSB of next vector
				xDataX <<= 1; // LSB = 0, discard MSB
			}
			index = ((quadrant<<1)&(~(0-signum)))+offset;
			system[index] ^= xDataX;
			index = ((quadrant<<1)&(0-signum))+offset;
			system[index] ^= thisDataX;
			system[index+quadrant] = thisDataY ^ thisDataX;
			system[yIndexY] = yDataY ^ thisDataX;
		}
	}
}

template<class Rng>
void Kpz::SchedulerCPUCBBit<Rng>::Func_MCSLocalBase::pqSyncRngUpdate(int yMin, int ySup)
{
	for(int y = yMin; y < ySup; ++y)
		for(int x = 0; x < m_this->m_size.dimX(); ++x)
		{
			if(((x^y)&1) != m_this->m_signum)
				continue;

			Slope cx = m_this->m_size.slopeX_CBBit(x, y);
			Slope cy = m_this->m_size.slopeY_CBBit(x, y);
			Slope r = m_this->m_size.slopeX_CBBit((x+1)&m_this->m_sizeCache.mDimX(), y);
			Slope l = m_this->m_size.slopeY_CBBit(x, (y+1)&m_this->m_sizeCache.mDimY());

			const int config = cx(m_this->m_stash) | (cy(m_this->m_stash)<<1)
				| (r(m_this->m_stash)<<2) | (l(m_this->m_stash)<<3);

			const float rnd = m_rng->randomFloat();
			if(config == 0b1100)
			{ // p
				if(rnd >= m_this->m_disorderP[0])
					continue;
			}
			else if(config == 0b0011)
			{ // q
				if(rnd >= m_this->m_disorderQ[0])
					continue;
			}
			else continue;
			cx.flip(m_this->m_stash);
			cy.flip(m_this->m_stash);
			r.flip(m_this->m_stash);
			l.flip(m_this->m_stash);
		}
}

template<class Rng>
void Kpz::SchedulerCPUCBBit<Rng>::Func_MCSLocalBase::disorder2SyncRngUpdate(int yMin, int ySup)
{
	for(int y = yMin; y < ySup; ++y)
		for(int x = 0; x < m_this->m_size.dimX(); ++x)
		{
			if(((x^y)&1) != m_this->m_signum)
				continue;

			Slope cx = m_this->m_size.slopeX_CBBit(x, y);
			Slope cy = m_this->m_size.slopeY_CBBit(x, y);
			Slope r = m_this->m_size.slopeX_CBBit((x+1)&m_this->m_sizeCache.mDimX(), y);
			Slope l = m_this->m_size.slopeY_CBBit(x, (y+1)&m_this->m_sizeCache.mDimY());

			const int config = cx(m_this->m_stash) | (cy(m_this->m_stash)<<1)
				| (r(m_this->m_stash)<<2) | (l(m_this->m_stash)<<3);

			const float rnd = m_rng->randomFloat();
			if(config == 0b1100)
			{ // p
				if(rnd >= m_this->m_disorderP[cx(m_this->m_disorder)])
					continue;
			}
			else if(config == 0b0011)
			{ // q
				if(rnd >= m_this->m_disorderQ[cx(m_this->m_disorder)])
					continue;
			}
			else continue;
			cx.flip(m_this->m_stash);
			cy.flip(m_this->m_stash);
			r.flip(m_this->m_stash);
			l.flip(m_this->m_stash);
		}
}
