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

#ifndef KMC_INCLUDE_MAP_H
#define KMC_INCLUDE_MAP_H

#include <iostream>

namespace Kmc
{

/*! \brief Template for an \p N dimensional array.
 *
 * \details Abstracts efficient accesses by coordinates.
 * \tparam N spacial dimension
 * \tparam Dim Type to use for lateral sizes. A type restricting sizes to
 * powers two in exchange for fast bit-operations can be used here. Type has to
 * implement `operator*`.
 * \tparam Elem type of elements
 * \ingroup DGutils
 */
template<int N, class Dim, class Elem = double>
class Map
{
	public:

	/*! \brief Coordinate vector.
	 * \details For `N == 3` compatible to ::Point.
	 */
	struct Point
	{
		int coord[N];

			Point() {}
		/*! \param v initial value for all coordinates. */
			Point (int v) {
			for(int a = 0; a < N; ++a)
				coord[a] = v;
		}

		/*! \param i component index, periodically mapped to 0...N-1
		 * \return respective component
		 */
		const int& operator[](int i) const {return coord[i%N];}
		/*! \param i component index, periodically mapped to 0...N-1
		 * \return reference to respective component
		 */
		int& operator[](int i) {return coord[i%N];}

		/*! \param map Map to access
		 * \return corresponding index of element in \p map
		 */
		size_t index(Map<N, Dim, Elem>& map) const {
			size_t i = coord[N-1];
			for(int a = N-2; a >= 0; --a)
			{
				i *= map.m_dim[a];
				i += coord[a];
			}
			return i;
		}

		/*! \return square of Cartesian norm */
		double abs2() const {
			double ret = coord[0]*(double)coord[0];
			for(int a = 1; a < N; ++a)
				ret += coord[a]*(double)coord[a];
			return ret;
		}
	};

	protected:
	Elem* m_data;
	Dim m_dim[N];
	size_t m_size, m_stride[N];

	/*! Allocator predicate, see setup(Alloc).
	 * \param size number of elements
	 * \return pointer to newly allocated memory to store the specified number
	 * of elements.
	 */
	static Elem* allocNew(size_t size) {return new Elem[size];}

	/*! Set up map with set parameters.
	 * \tparam Alloc Predicate to allocate the internal data array.
	 * \param alloc instance of \p Alloc
	 */
	template<class Alloc = Elem* (size_t)>
	void setup(Alloc allocator = allocNew) {
		m_stride[0] = 1;
		m_size = (size_t)m_dim[0];
		for(int a = 1; a < N; ++a)
		{
			m_stride[a] = m_stride[a-1] * m_dim[a-1];
			m_size *= m_dim[a];
		}
		m_data = allocator(m_size);
	}

	public:

	/*! Construct cubic map with given lateral size.
	 * \param dim lateral size
	 */
		Map(const Dim& dim) {
			for(int a = 0; a < N; ++a)
				m_dim[a] = dim;
			setup();
		}
	/*! Construct map with given dimensions.
	 * \param dim dimensions
	 */
		Map(const Point& dim) {
			for(int a = 0; a < N; ++a)
				m_dim[a] = Dim(dim.coord[a]);
			setup();
		}
	/*! Construct map with given dimensions.
	 * \tparam Alloc Predicate to allocate the internal data array.
	 * \param dim dimensions (here array of Dim)
	 * \param alloc instance of \p Alloc
	 */
		template<class Alloc = Elem* (size_t)>
		Map(const Dim dim[N], Alloc allocator = allocNew) {
			for(int a = 0; a < N; ++a)
				m_dim[a] = dim[a];
			setup(allocator);
		}
		~Map() {
			delete[] m_data;
		}

	/*! Access element by coordinates. */
	Elem& operator[] (const Point& i) {return m_data[i.index(*this)];}
	/*! Access element by index. */
	Elem& operator[] (const int index) {return m_data[index];}
	/*! \overload operator(const Point&) */
	const Elem& operator[] (const Point& i) const {return m_data[i.index(*this)];}
	/*! \overload operator(const int) */
	const Elem& operator[] (const int index) const {return m_data[index];}
	/*! Access element by coordinates. */
	Elem& at (const Point& i) {return m_data[i.index(*this)];}
	/*! Access element by index. */
	Elem& at (const int index) {return m_data[index];}
	/*! \overload at(const Point&) */
	const Elem& at (const Point& i) const {return m_data[i.index(*this)];}
	/*! \overload at(const int) */
	const Elem& at (const int index) const {return m_data[index];}

	/*! \return number of elements */
	size_t size() const {return m_size;}
	/*! \param n dimension index
	 * \return element stride in given dimension.
	 */
	size_t stride(int n) const {return m_stride[n];}
	/*! \param n dimension index
	 * \return size in given dimension.
	 */
	Dim dim(int n) const {return m_dim[n];}
	const Dim* dim() const {return m_dim;}
	/*! \return pointer to internal array */
	Elem* data() {return m_data;}
	/*! Disown internal array, do not delete.
	 * \warning The internal data pointer is set to null by this operation,
	 * thus subsequent accesses will be invalid. Only use immediately before
	 * deleting the map if the data is to be preserved.
	 */
	void dropData() {m_data = 0;}
};

/*! Predicate that writes a value to the given stream.
 * \tparam type of value
 * \param o ostream
 * \param t value to write
 * \return \p o
 * \relates Kmc::Map
 */
template<class T>
inline std::ostream& outOperator(std::ostream& o, const T& t) {return o << t;}

/*! Write contents of a 2d  map to stream in gnuplots pm3d format.
 * \tparam Dim Dim parameter of map (does not have to be specified explicitly)
 * \tparam Elem element type of map (does not have to be specified explicitly)
 * \tparam Output print predicate
 * \param map 2d map to write
 * \param o ostream to write to
 * \param cellDim Lateral size to assume of each element, i.e. scale
 * coordinates by this factor.
 * \param out instance of Output
 * \relates Kmc::Map
 */
template<class Dim, class Elem, class Output = std::ostream& (std::ostream&, const Elem&)>
void writePM3D(const Map<2, Dim, Elem> &map, std::ostream& o, double cellDim = 1., Output out = outOperator<Elem>);
/*! Write contents of a map to stream as ASCII table
 * \tparam N spacial dimension map (does not have to be specified explicitly)
 * \tparam Dim Dim parameter of map (does not have to be specified explicitly)
 * \tparam Elem element type of map (does not have to be specified explicitly)
 * \tparam Output print predicate
 * \param map 2d map to write
 * \param o ostream to write to
 * \param cellDim Lateral size to assume of each element, i.e. scale
 * coordinates by this factor.
 * \param out instance of Output
 */
template<int N, class Dim, class Elem, class Output = std::ostream& (std::ostream&, const Elem&)>
void writeASCII(const Map<N, Dim, Elem> &map, std::ostream& o, double cellDim = 1., Output out = outOperator<Elem>);
/*! Write contents of a 2d  map to stream in as ASCII matrix readable by
 * GnuOctave
 * \tparam Dim Dim parameter of map (does not have to be specified explicitly)
 * \tparam Elem element type of map (does not have to be specified explicitly)
 * \tparam Output print predicate
 * \param map 2d map to write
 * \param o ostream to write to
 * \param name name the matrix will have when loaded with Octave
 * \param out instance of Output
 * \relates Kmc::Map
 */
template<class Dim, class Elem, class Output = std::ostream& (std::ostream&, const Elem&)>
void writeOctaveMat(const Map<2, Dim, Elem> &map, std::ostream& o, const char* name = "a", Output out = outOperator<Elem>);
/*! Compute variance of map's elements.
 * \tparam N spacial dimension map (does not have to be specified explicitly)
 * \tparam Dim Dim parameter of map (does not have to be specified explicitly)
 * \tparam Elem element type of map (does not have to be specified explicitly)\n
 * Required to implement `operator+=`, `operator*` and `operator/(int or
 * float)` as well as a constructor from int or float. (... or double)
 * \param map Map
 * \return variance of given map's elements.
 * \relates Kmc::Map
 */
template<int N, class Dim, class Elem>
Elem variance(const Map<N, Dim, Elem> &map);
}

#include "mapTemplate.cpp"
#endif
