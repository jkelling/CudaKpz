/***************************************************************************
*   Copyright 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

template<class T>
template<class Iter>
void HistogramMoments<T>::set(Iter begin, Iter end)
{
	m_samples = 0;
	m_mean = m_mean2 = 0;
	int n = 0;
	for(auto a = begin; a != end; ++a, ++n)
	{
		m_samples += *a;
		const auto tmp = (*a)*n;
		m_mean += tmp;
		m_mean2 += tmp*n;
	}
	m_mean /= m_samples;
	m_mean2 /= m_samples;

	{
		m_median = 0;
		int tmpSum = 0;
		for(auto a = begin; a != end; ++a, ++m_median)
		{
			tmpSum += *a;
			if(tmpSum*2 >= m_samples)
				break;
		}
	}
}
