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

#include "HistogramMoments.h"
#include "binnedVector.h"

template<class T>
class BinnedVectorMoments : public HistogramMoments<T>
{
	T m_bin;

	public:
		BinnedVectorMoments(const BinnedVector<T>& vec)
			: HistogramMoments<T>(vec.begin(), vec.end()), m_bin(vec.bin()) {}

	void set(const BinnedVector<T>& vec) {
		m_bin=vec.bin;
		set(vec.begin(), vec.end());
	}
	void setBin(T bin) {m_bin=bin;}

	T median() const {return HistogramMoments<T>::median()*m_bin;}
	T mean() const {return HistogramMoments<T>::mean()*m_bin;}
	T mean2() const {return HistogramMoments<T>::mean2()*m_bin*m_bin;}
	T bin() const {return bin;}
};
