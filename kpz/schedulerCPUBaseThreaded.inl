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

inline void Kpz::SchedulerCPUBaseThreaded::pqOneUpdate(int x, int y)
{
	LatticePoint c = m_size.latticePoint(x, y);
	Slope r = m_size.slopeX((x+1)&m_sizeCache.mDimX(), y);
	Slope l = m_size.slopeY(x, (y+1)&m_sizeCache.mDimY());
	
	if((!c(m_stash)) && r(m_stash) && l(m_stash))
	{
		c.flip(m_stash);
		r.flip(m_stash);
		l.flip(m_stash);
	}
}
