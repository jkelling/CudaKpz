/***************************************************************************
*   Copyright 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

template<class T>
class HistogramMoments
{
	T m_median, m_samples;
	T m_mean, m_mean2;

	public:
	template<class Iter>
		HistogramMoments(Iter begin, Iter end) {set(begin,end);}

	template<class Iter>
	void set(Iter begin, Iter end);

	T median() const {return m_median;}
	T mean() const {return m_mean;}
	T samples() const {return m_samples;}
};

#include "HistogramMoments.inl"
