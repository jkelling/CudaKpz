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


template<int N, class Dim>
void HeightMap<N, Dim>::tryNeighbours(const size_t index, const System<N, Dim>& system, const bool min) {
	LatticePoint<N, Dim> l(dim(), index); 
#ifdef KPZ_DBG_LOG
	std::cerr << l << " extr " << m_system[index] << '\n';
#endif
	const float neighbourVal = m_system[index] + (min ? 1 : -1);
	for(int a = N-1; a >= 0; --a)
	{
		size_t t = incCoord(l, a, index);
		if(isnan(m_system[t]))
		{
			m_system[t] = neighbourVal;
			const int him = system.get(t);
#ifdef KPZ_DBG_LOG
			std::cerr << l << " extr -> " << LatticePoint<N>(dim(), t) << ' ' << him << '\n'; 
#endif
			if(him == 0) {
				tryNeighbours(t, system, true);
			}
			else if(him == 2*N) {
				tryNeighbours(t, system, false);
			}
			else {
				tryNeighbours(t, system);
			}
		}
		t = decCoord(l, a, index);
		if(isnan(m_system[t]))
		{
			m_system[t] = neighbourVal;
			const int him = system.get(t);
#ifdef KPZ_DBG_LOG
			std::cerr << l << " extr -> " << LatticePoint<N>(dim(), t) << ' ' << him << '\n'; 
#endif
			if(him == 0) {
				tryNeighbours(t, system, true);
			}
			else if(him == 2*N) {
				tryNeighbours(t, system, false);
			}
			else {
				tryNeighbours(t, system);
			}
		}
	}
}

template<int N, class Dim>
void HeightMap<N, Dim>::tryNeighbours(const size_t index, const System<N, Dim>& system) {
	LatticePoint<N, Dim> l(dim(), index); 
#ifdef KPZ_DBG_LOG
	std::cerr << l << " normal " << m_system[index] << '\n';
#endif
	{
		int nLow = 0, nHigh = 0;
		{
			const float val= m_system[index];
			for(int a = N-1; a >= 0; --a)
			{
				size_t t = incCoord(l, a, index);
				if(isnan(m_system[t]))
				{
					const int him = system.get(t);
					if(him == 0) {
						++nLow;
						m_system[t] = val- 1;
						tryNeighbours(t, system, true);
					}
					else if(him == 2*N) {
						++nHigh;
						m_system[t] = val+ 1;
						tryNeighbours(t, system, false);
					}
				}
				else
				{
					if(m_system[t] > val)
						++nHigh;
					else
						++nLow;
				}

				t = decCoord(l, a, index);
				if(isnan(m_system[t]))
				{
					const int him = system.get(t);
					if(him == 0) {
						++nLow;
						m_system[t] = val- 1;
						tryNeighbours(t, system, true);
					}
					else if(him == 2*N) {
						++nHigh;
						m_system[t] = val+ 1;
						tryNeighbours(t, system, false);
					}
				}
				else
				{
					if(m_system[t] > val)
						++nHigh;
					else
						++nLow;
				}
			}
		}
		if(nHigh + nLow == 2*N)
			return;
		const int me = system.get(index);
#ifdef KPZ_DBG_LOG
		std::cerr << " normal nLow " << nLow << " nHigh " << nHigh << " me " << me << '\n';
#endif
		if(nLow + 2*N - me == 2*N) {
			tryNeighbours(index, system, true);
			return;
		}
		else if(nHigh + me == 2*N) {
			tryNeighbours(index, system, false);
			return;
		}
	}

	for(int a = N-1; a >= 0; --a)
	{
		size_t t = incCoord(l, a, index);
		if(isnan(m_system[t])) {
			trySet(t, system);
		}
		t = decCoord(l, a, index);
		if(isnan(m_system[t])) {
			trySet(t, system);
		}
	}
}

template<int N, class Dim>
void HeightMap<N, Dim>::trySet(const size_t index, const System<N, Dim>& system) {
	LatticePoint<N, Dim> l(dim(), index); 
#ifdef KPZ_DBG_LOG
	std::cerr << l << " set " << m_system[index] << '\n';
#endif
	const int me = system.get(index);
	float neighbourVal = NAN; // save lower neighbour value in neighbourVal
	int nLow = 0;
	for(int a = N-1; a >= 0; --a)
	{
		float t = m_system[incCoord(l, a, index)];
		if(isnan(t));
		else if(isnan(neighbourVal))
		{
#ifdef KPZ_DBG_LOG
			std::cerr << " set initial nv " << neighbourVal << " t " << t << '\n';
#endif
			neighbourVal = t;
			nLow = 1;
		}
		else if(t != neighbourVal)
		{
#ifdef KPZ_DBG_LOG
			std::cerr << " set final t " << t << " nv " << neighbourVal << '\n'; 
#endif
			if(t < neighbourVal)
				--neighbourVal;
			else
				++neighbourVal;
			nLow = -1;
			break;
		}
		else
			++nLow;

		t = m_system[decCoord(l, a, index)];
#ifdef KPZ_DBG_LOG
		std::cerr << " set param t " << t << " index " << incCoord(l, a, index) << '\n';
#endif
		if(isnan(t));
		else if(isnan(neighbourVal))
		{
#ifdef KPZ_DBG_LOG
			std::cerr << " set initial nv " << neighbourVal << " t " << t << '\n';
#endif
			neighbourVal = t;
			nLow = 1;
		}
		else if(t != neighbourVal)
		{
#ifdef KPZ_DBG_LOG
			std::cerr << " set final  t " << t << " nv " << neighbourVal << '\n'; 
#endif
			if(t < neighbourVal)
				--neighbourVal;
			else
				++neighbourVal;
			nLow = -1;
			break;
		}
		else
			++nLow;
	}
	if(nLow == -1);
	else if(me + nLow > (2*N))
		++neighbourVal;
	else if(me - nLow < 0)
		--neighbourVal;
	else {
#ifdef KPZ_DBG_LOG
		std::cerr << " set fail\n";
#endif
		return;
	}
#ifdef KPZ_DBG_LOG
	std::cerr << " set " << nLow << ' ' << neighbourVal << '\n';
#endif
	m_system[index] = neighbourVal;
	tryNeighbours(index, system);
}
