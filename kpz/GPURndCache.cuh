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

#pragma once

/*! \brief produces random bits insteat of float, p=0.5
 *
 * Should be used in synchronized fashion by threads of a warp to eliminate
 * negative effect of branches.
 */
template<class Rng>
class BinaryRndFloatCache_5 : public Rng
{
	public:

		using Rng::Rng;

	class Device
	{
		protected:
		typename Rng::Device rng;
		int m_cache;
		unsigned char m_cacheFill;

		__device__ void fillCache() {
			m_cache = rng.random();
		}

		public:

		inline __device__ int random() {return rng.random();}
		inline __device__ float randomFloat() {
			if(m_cacheFill)
			{
				m_cache>>=1;
				--m_cacheFill;
			}
			else
			{
				fillCache();
				m_cacheFill = 31; //! \todo we could ask the rng about the number of random bits
			}
			return m_cache&1;
		}

		__device__ Device(void* d) : rng(d), m_cacheFill(0) {}
	};
};
