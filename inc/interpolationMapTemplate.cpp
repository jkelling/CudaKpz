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

#include <cmath>

/*!
 * \brief Internal helper for InterpolationMap.
 * \details Holds mediation parameters.
 */
template<int N, class Dim, class T>
struct MediationParam
{
	typename Kmc::InterpolationMap<N,Dim,T>::PointF up, down;
	double weight;

	MediationParam() : up(1.), down(1.), weight(1.) {}
};

template<int N, class Dim, class T>
template<class Mediator>
bool Kmc::InterpolationMap<N, Dim, T>::add(const T& v, PointF p, Mediator mediator)
{
	Point i;
	MediationParam<N,Dim,T> mediation;
	int includeMask = 0; //bitmask for directions that must be included when computing share
	for(int a = 0; a < N; ++a)
	{
		p.coord[a] /= m_cellDim;
		i.coord[a] = round(p.coord[a]);
		if(i.coord[a] < -1 || i.coord[a] > (int)m_dim[a])
				return false;
		else if(i.coord[a] < 0 || i.coord[a] >= (int)m_dim[a])
		{
			includeMask |= (1<<a);
		}
		p.coord[a] -= i.coord[a]-.5;
		const double b = p.coord[a] / m_cellBorder;
		if(b < 1.)
		{
			if(i.coord[a] > 0)
				mediation.down.coord[a] = mediator(b);
		}
		else if(b > m_cellBorderHigh && (i.coord[a] < ((int)m_dim[a]-1)))
			mediation.up.coord[a] =  mediator(1. - (b - m_cellBorderHigh));
	}

	const size_t index = i.index(*this);
	for(int n = (1<<N)-1; n >= includeMask; --n)
	{
		if((n&includeMask) != includeMask)
			continue; //this neighbor does not lie within the system
		size_t t = index;
		mediation.weight = 1.;
		for(int a = 0; a < N; ++a)
		{
			if(mediation.down.coord[a] < 1.)
			{
				if((n>>a)&1)
				{
					mediation.weight *= (1. - mediation.down.coord[a]);
					t -= m_stride[a];
				}
				else
					mediation.weight *= (mediation.down.coord[a]);
			}
			else if(mediation.up.coord[a] < 1.)
			{
				if((n>>a)&1)
				{
					mediation.weight *= (1. - mediation.up.coord[a]);
					t += m_stride[a];
				}
				else
					mediation.weight *= (mediation.up.coord[a]);
			}
			else if((n>>a)&1)
			{
				goto noShare; // this cell does not have a share
			}
		}
		m_data[t].add(v, mediation.weight);
noShare:;
	}
	return true;
}

template<int N, class Dim, class T>
void Kmc::InterpolationMap<N, Dim, T>::finish()
{
	for(size_t a = 0; a < m_size; ++a)
		if(m_data[a].weights > 0.)
			m_data[a].weight(1.);
}

/*!
 * \brief Internal helper for InterpolationMap.
 * \details Part of TMP for InterpolationMap::interpolateEmpty().
 */
template<int n, int N, class Dim, class T, class Interpolator>
struct InterpolateEmptyHelper
{
	static void iterate(Kmc::Map<N,Dim,Kmc::InterpolationItem<T> >& map
		, typename Kmc::Map<N,Dim,Kmc::InterpolationItem<T> >::Point& pos, size_t& index, Interpolator interpolator)
	{
		for(pos.coord[n] = 0; pos.coord[n] < (int)map.dim(n); ++pos.coord[n])
		{
			InterpolateEmptyHelper<n-1, N, Dim, T, Interpolator>::iterate(map, pos, index, interpolator);
		}
	}
};

/*!
 * \brief Internal helper for InterpolationMap.
 * \details Part of TMP for InterpolationMap::interpolateEmpty().
 */
template<int N, class Dim, class T, class Interpolator>
struct InterpolateEmptyHelper<0, N, Dim, T, Interpolator>
{
	static void iterate(Kmc::Map<N,Dim,Kmc::InterpolationItem<T> >& map
		, typename Kmc::Map<N,Dim,Kmc::InterpolationItem<T> >::Point& pos, size_t& index, Interpolator interpolator)
	{
		for(pos.coord[0] = 0; pos.coord[0] < (int)map.dim(0); ++pos.coord[0], ++index)
		{
			if(map[index].weights > 0.)
				continue;
			for(int s = 0; s < N; ++s)
			{
				const bool infoLower = (pos.coord[s] == 0) ? false : map[index-map.stride(s)].weights > 0.;
				const bool infoUpper = (pos.coord[s] >= ((int)map.dim(s)-1)) ? false : map[index+map.stride(s)].weights > 0.;
				if((!infoLower) & infoUpper)
					map[index] += map[index+map.stride(s)];
				else if(infoLower & (!infoUpper))
					map[index] += map[index-map.stride(s)];
				else if(infoLower & infoUpper)
				{
					map[index] += interpolator(map[index-map.stride(s)], map[index+map.stride(s)]);
				}
			}
			if(map[index].weights > 0.)
				map[index].weight(0.);
		}
	}
};

template<int N, class Dim, class T>
template<class Interpolator>
void Kmc::InterpolationMap<N, Dim, T>::interpolateEmpty(Interpolator interpolator)
{
	size_t index = 0;
	Point pos;
	InterpolateEmptyHelper<N-1, N, Dim, T, Interpolator>::iterate(*this, pos, index, interpolator);
}
