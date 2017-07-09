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

#ifndef KPZ_TIMER_H
#define KPZ_TIMER_H

#include "ctime"

namespace Kmc
{
	class Timer
	{
		time_t m_end;

		public:

			Timer(time_t t) : m_end(t) {}
			Timer(float h) {set(h);}

		void set(float th);
		void set(time_t t) {m_end = t;}

		time_t end() const {return m_end;};
		bool atEnd() const {return m_end <= time(0);}
	};

	class TimerSingleton
	{
		static Timer* m_instance;

			TimerSingleton();
			TimerSingleton(const TimerSingleton&);
		public:
		static void create(time_t t) {
			if(m_instance)
				m_instance->set(t);
			else
				m_instance = new Timer(t);
		}

		static void create(float th) {
			if(m_instance)
				m_instance->set(th);
			else
				m_instance = new Timer(th);
		}

		static void destroy() {
			delete m_instance;
			m_instance = 0;
		}

		static bool active() {
			return m_instance != 0;
		}

		static bool atEnd() {
			return m_instance != 0 && m_instance->atEnd();
		}
	};
}

#endif
