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

#ifndef N_KPZ_DIM_H
#define N_KPZ_DIM_H

#include "../../kmc/refactorized/lattice.h"

namespace nKPZ
{

class IntDim
{
	protected:
	int m_val;

	public:

	inline IntDim() : m_val(0) {}
	inline IntDim(int val) : m_val(val) {}
	
	inline IntDim& operator= (int val) {m_val = val; return *this;}
	inline IntDim& operator= (size_t val) {m_val = val; return *this;}
	inline IntDim& operator= (const IntDim& val) {m_val = val.m_val; return *this;}
	inline int operator* (int val) const {return m_val * val;}
	inline size_t operator* (size_t val) const {return m_val * val;}
	inline double operator* (double val) const {return m_val * val;}

	explicit inline operator int () const {return m_val;}
	explicit inline operator size_t () const {return m_val;}

	class Mod {
		protected:
		int m_val;

		public:
			
		inline Mod(const IntDim& dim) : m_val(dim.m_val) {}

		inline int operator() (int a) const {return a%m_val;}
		inline size_t operator() (size_t a) const {return a%m_val;}
	};

	inline bool isEven() const {return ~m_val&1;}

};

inline int operator/ (int a, const IntDim& b) {return a / (int)b;}
inline size_t operator/ (size_t a, const IntDim& b) {return a / (int)b;}
inline int operator/= (int& a, const IntDim& b) {return a /= (int)b;}
inline size_t operator/= (size_t& a, const IntDim& b) {return a /= (int)b;}
inline int operator% (int a, const IntDim& b) {return a % (int)b;}
inline size_t operator% (size_t a, const IntDim& b) {return a % (int)b;}
inline int operator* (int a, const IntDim& b) {return b*a;}
inline size_t operator* (size_t a, const IntDim& b) {return b*a;}
inline double operator* (double a, const IntDim& b) {return b*a;}
inline int operator*= (int& a, const IntDim& b) {return a=b*a;}
inline size_t operator*= (size_t& a, const IntDim& b) {return a=b*a;}

class LogDim
{
	protected:
	int m_val;

	public:

	inline LogDim() : m_val(0) {}
	inline LogDim(int val) : m_val(Kmc::log2up(val)) {}
	
	inline LogDim& operator= (int val) {m_val = val; return *this;}
	inline LogDim& operator= (size_t val) {m_val = val; return *this;}
	inline LogDim& operator= (const LogDim& val) {m_val = val.m_val; return *this;}
	inline int operator* (int val) const {return  val << m_val;}
	inline size_t operator* (size_t val) const {return val << m_val;}
	inline double operator* (double val) const {return (1<<m_val) * val;}

	explicit inline operator int () const {return (1<<m_val);}
	explicit inline operator size_t () const {return (1<<m_val);}

	class Mod {
		protected:
		int m_val;

		public:
			
		inline Mod(const LogDim& dim) : m_val((1<<dim.m_val)-1) {}

		inline int operator() (int a) const {return a&m_val;}
		inline size_t operator() (size_t a) const {return a&m_val;}
	};

	inline bool isEven() const {return m_val;}

	friend inline int operator/ (int a, const LogDim& b);
	friend inline size_t operator/ (size_t a, const LogDim& b);
	friend inline int operator/= (int& a, const LogDim& b);
	friend inline size_t operator/= (size_t& a, const LogDim& b);
	friend inline int operator% (int a, const LogDim& b);
	friend inline size_t operator% (size_t a, const LogDim& b);

};

inline int operator/ (int a, const LogDim& b) {return a >> b.m_val;}
inline size_t operator/ (size_t a, const LogDim& b) {return a >> b.m_val;}
inline int operator/= (int& a, const LogDim& b) {return a >>= b.m_val;}
inline size_t operator/= (size_t& a, const LogDim& b) {return a >>= b.m_val;}
inline int operator% (int a, const LogDim& b) {return a &((1<<b.m_val)-1);}
inline size_t operator% (size_t a, const LogDim& b) {return a &((1<<b.m_val)-1);}
inline int operator* (int a, const LogDim& b) {return b*a;}
inline size_t operator* (size_t a, const LogDim& b) {return b*a;}
inline double operator* (double a, const LogDim& b) {return b*a;}
inline int operator*= (int& a, const LogDim& b) {return a=b*a;}
inline size_t operator*= (size_t& a, const LogDim& b) {return a=b*a;}

}

#endif
