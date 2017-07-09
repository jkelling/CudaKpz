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

#ifndef KMC_INCLUDE_INTERPOLATION_MAP
#define KMC_INCLUDE_INTERPOLATION_MAP

#include "map.h"

namespace Kmc
{

/*! \brief Element type for InterpolationMap.
 *
 * \details Wraps desired element type for use with InterpolationMap. Hold
 * weight of the element.
 * \tparam T element type
 */
template<class T>
struct InterpolationItem
{
	T value;
	double weights;

		InterpolationItem() : value(), weights(0.) {}
	/*! \param v element value
	 * \param w element weight or sum of weights
	 */
		InterpolationItem(const T& v, double w) : value(v), weights(w) {}

	/*! Add another weighted value \p v to this element.
	 * \param v value
	 * \param weight weight of the added value
	 */
	inline void add(const T& v, double weight) {
		weights += weight;
		value += v*weight;
	}
	
	/*! Add two InterpolationItem\e s. Add weights and values
	 * \param i other InterpolationItem
	 * \return reference to *this
	 */
	inline InterpolationItem& operator+= (const InterpolationItem& i) {
		weights += i.weights;
		value += i.value;
		return *this;
	}

	/*! Set new weight of the element. Apply accumulated weights to value
	 * beforehand (weighted average).
	 * \param weight new weight
	 */
	inline void weight(double weight) {
		value /= weights;
		weights = weight;
	}

	/*! Predicate that writes the elements value to the given stream.
 	 * \param o ostream
	 * \param i InterpolationItem
	 * \return \p o
	 */
	static inline std::ostream& printValue(std::ostream& o, const InterpolationItem& i) {
		return o << i.value;
	}
};

/*! \brief Template for an \p N dimensional interpolated map.
 *
 * \details Allows writes to floating point coordinates distributing the
 * written value among neighboring elements.
 * \tparam N spacial dimension
 * \tparam Dim Type to use for lateral sizes. A type restricting sizes to
 * powers two in exchange for fast bit-operations can be used here. Type has to
 * implement `operator*`.
 * \tparam Elem type of elements
 * \ingroup DGutils
 */
template<int N, class Dim, class T = double>
class InterpolationMap : public Map<N, Dim, InterpolationItem<T> >
{
	using Map<N, Dim, InterpolationItem<T> >::m_size;
	using Map<N, Dim, InterpolationItem<T> >::m_dim;
	using Map<N, Dim, InterpolationItem<T> >::m_stride;
	using Map<N, Dim, InterpolationItem<T> >::m_data;

	//! Lateral cell size.
	double m_cellDim;
 	//! Fraction of lateral cell size that is considered shared.
	double m_cellBorder;
	/*! There is a trick connected to the value assigned to this member: Assume
	 * `f = x / m_cellDim`, where \c x is a component of a point. Then if `f <
	 * m_cellBorder` the point is shared with the lower neighbor, if `f >
	 * m_cellBorderHigh` it is shared with the upper neighbor. This is why
	 * `m_cellBorderHigh = 1./cellBorder-1.`
	 */
	double m_cellBorderHigh;

	public:

	/*! Mediation predicate: Always share equally (.5 for both cells).
	 * \param r relative position of point on shared border
	 * \return .5
	 */
	static double MediatorConst(const double r) {return .5;}
	/*! Mediation predicate: Mediate linearly.
	 * \param r relative position of point on shared border
	 * \return linear interpolation `r*.5 + .5`\n (`r==0` means the point exactly
	 * in the middle between both cells)
	 */
	static double MediatorLinear(const double r) {return r*.5 + .5;}

	/*! Interpolation predicate: Linear interpolation.
	 * \param a neighbor
	 * \param b facing neighbor
	 * \return average of both value and weight
	 */
	static InterpolationItem<T> InterpolatorLinear(const InterpolationItem<T>& a, const InterpolationItem<T>& b) {
		return InterpolationItem<T>((a.value+b.value)/2., (a.weights+b.weights)/2.);
	}

	typedef typename Map<N, Dim, InterpolationItem<T> >::Point Point;

	/*! \brief Floating point coordinate vector.
	 * \details For `N == 3` compatible to ::PointD.
	 */
	struct PointF
	{
		double coord[N];

			PointF() {}
		/*! \param v initial value for all coordinates. */
			PointF (double v) {
			for(int a = 0; a < N; ++a)
				coord[a] = v;
		}
	};

	/*! Construct cubic map.
	 * \param dim lateral size
	 * \param cellDim lateral size to assign to cells
	 * \param cellBorder Fraction of cell to be considered shared with neighboring cells (in (0,.5)).
	 */
		InterpolationMap(Dim dim, double cellDim, double cellBorder)
			: Map<N, Dim, InterpolationItem<T> >(dim)
			, m_cellDim(cellDim), m_cellBorder(cellBorder), m_cellBorderHigh(1./cellBorder-1.) {}
	/*! Construct cubiod map.
	 * \param dim dimensions
	 * \param cellDim lateral size to assign to cells
	 * \param cellBorder Fraction of cell to be considered shared with neighboring cells (in (0,.5)).
	 */
		InterpolationMap(const Dim dim[N], double cellDim, double cellBorder)
			: Map<N, Dim, InterpolationItem<T> >(dim)
			, m_cellDim(cellDim), m_cellBorder(cellBorder), m_cellBorderHigh(1./cellBorder-1.) {}

	/*! Add value \p v to cells sharing given coordinates, mediated using the
	 * predicate \p Mediator.
	 * \tparam Mediator Predicate for computing the weights of the value shall
	 * have in two sharing cells. The predicate takes a value \c r in (0,1)
	 * stating how far the coordinate are inside the primary cell, 0 meaning
	 * the coordinates lie exactly at the cell boundary. The return value \c w
	 * is the weight of \p v for the primary cell, \c w - 1 being the weight
	 * for the other cell.
	 * \param v value
	 * \param p float coordinates
	 * \param mediator instance of Mediator
	 * \return whether \p is in range, thus if the value was added
	 */
	template<class Mediator = double(double)>
	bool add(const T& v, PointF p, Mediator mediator = MediatorLinear);
	/*! Compute weighted average for all cells: Call
	 * InterpolationItem::weight(double) (with 1.) for all cell that are not
	 * empty (i.e.  `weight != 0`).
	 */
	void finish();
	/*! Interpolate values for empty cells (those with `weight == 0`) from
	 * non-empty neighbors. Cells that do not have non-empty neighbors stay
	 * empty.
	 *
	 * The interpolation is done for each axis separately, then
	 * weighted average of these results is taken. Form axes where there are no
	 * non-empty neighbors there is no contribution.
	 * \tparam Interpolator Predicate interpolation between two facing neighbors.
	 * \param interpolator instance of Interpolator
	 */
	template<class Interpolator = double(const InterpolationItem<T>&, const InterpolationItem<T>&)>
	void interpolateEmpty(Interpolator interpolator = InterpolatorLinear); // does not use interpolated data for interpolation
	//! \todo Out of place version of interpolateEmpty().

	/*! \return defined lateral cell size. */
	double cellDim() const {return m_cellDim;}
	/*! \return fraction of lateral cell size considered shared with neighbors.*/
	double cellBorder() const {return m_cellBorder;}
};

/*! \param o ostream to write ot
 * \param i InterpolationItem to write
 * \return \p o
 * \relates InterpolationItem
 */
template<class T>
std::ostream& operator<<(std::ostream& o, const Kmc::InterpolationItem<T>& i)
{
	return o << i.value << '\t' << i.weights;
}

/*! Write contents of a 2d InterpolationItem to stream in as ASCII matrix
 * readable by GnuOctave
 * \tparam Dim Dim parameter of map (does not have to be specified explicitly)
 * \tparam Elem element type of map (does not have to be specified explicitly)
 * \tparam Output print predicate
 * \param map 2d map to write
 * \param o ostream to write to
 * \param name name the matrix will have when loaded with Octave
 * \param out instance of Output
 * \relates Kmc::InterpolationMap
 */
template<class Dim, class Elem>
inline void writeOctaveMat(const InterpolationMap<2, Dim, Elem> &map, std::ostream& o, const char* name = "a")
{
	writeOctaveMat(map, o, name, InterpolationItem<Elem>::printValue);
}

}

#include "interpolationMapTemplate.cpp"
#endif
