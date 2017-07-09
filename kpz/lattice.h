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

#ifndef KPZ_LATTICE_H
#define KPZ_LATTICE_H

namespace Kpz
{
	struct LatticePoint;

	struct Slope
	{ // manipulates single slope (one of the two bits of each site)
		int index, shift;

		inline Slope (int i, int s) 
			: index(i), shift(s)
		{}
		inline Slope () {}

		inline int operator() (const unsigned int* system) const {
			return (system[index] >> shift)&1;
		}
		
		inline void set (unsigned int* system) const {
			system[index] |= 1<<shift; //set slope
		}

		inline void unset (unsigned int* system) const {
			system[index] &= ~(1<<shift); //unset slope
		}

		inline void flip (unsigned int* system) const {
			system[index] ^= 1<<shift; //unset slope
		}

		inline void set (int v, unsigned int* system) const {
			// v must contain Bit0, rest nulled
			system[index] &= ~(1<<shift); //null slope
			system[index] |= v<<shift; //set value
		}

	};

	struct LatticePoint : public Slope
	{ // manipulates whole site (both slopes), assuming (shift&1==0)

		inline LatticePoint (const Slope& s)
			: Slope(s.index, s.shift&(~1))
		{}
		inline LatticePoint () {}
		inline LatticePoint (int i, int s) : Slope(i,s) {}

		inline int operator() (const unsigned int* system) const {
			return (system[index] >> shift)&3;
		}
		inline int operator() (const unsigned int system) const {
			return (system >> shift)&3;
		}
		
		inline void set (unsigned int* system) const {
			system[index] |= 3<<shift; //set site
		}

		inline void unset (unsigned int* system) const {
			system[index] &= ~(3<<shift); //unset site
		}

		inline void flip (unsigned int* system) const {
			system[index] ^= 3<<shift; //unset slope
		}

		inline void set (int v, unsigned int* system) const {
			// v must contain two bits, rest nulled
			system[index] &= ~(3<<shift); //null lattice site
			system[index] |= v<<shift; //set value
		}

	};

	//std::ostream& operator<< (std::ostream& o, const Slope& s) {
	//	return o << '(' << s.index << ", " << s.shift << ')';
	//}
}

#endif
