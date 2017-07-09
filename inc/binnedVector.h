/***************************************************************************
*   Copyright 2011 - 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#ifndef KMC_BINNED_VECTOR_H
#define KMC_BINNED_VECTOR_H

#include <vector>
#include <iostream>
#include <cmath>
#include <stdexcept>

template<class T>
class BinnedVector : public std::vector<T>
{
	double m_bin;

	public:

		BinnedVector(const double& bin = 1) : m_bin(bin) {}

	inline double bin() const {return m_bin;}
	inline void setBin(double bin) {m_bin = bin;}
	inline double binnedSize() const {return (this->size()-1)*m_bin;}

	inline T& operator[] (double i) {return std::vector<T>::operator[]((int)(i/m_bin));}
	inline const T& operator[] (double i) const {return std::vector<T>::operator[]((int)(i/m_bin));}

	inline void resize(double size) {std::vector<T>::resize((int)ceil(size/m_bin));}
	inline void resize(double size, const T& value) {std::vector<T>::resize((int)ceil(size/m_bin), value);}

	inline std::vector<T>& operator() () {return *this;}
	inline const std::vector<T>& operator() () const {return *this;}

};

template<class T>
class MultiColumnBinnedVector
{
	std::vector<BinnedVector<T> > m_y;
	double m_bin;

	public:

		MultiColumnBinnedVector(double bin = 1, int cols = 2) : m_bin(bin), m_y(0) {
		addCol(cols);
	}

	inline void clear() {m_y.clear();}

	inline double bin() const {return m_bin;}
	inline void setBin(double bin) {
		m_bin=bin;
		for(auto& col : m_y)
			col.setBin(m_bin);
	}

	inline std::vector<BinnedVector<T> >& list() {return m_y;}
	inline BinnedVector<T>& addCol() {
		m_y.emplace_back(m_bin);
		return m_y.back();
	}
	int addCol(int cols) {
		const int ret = m_y.size();
		m_y.reserve(ret+cols);
		for(int a = 0; a < cols; ++a)
			m_y.emplace_back(m_bin);
		return ret;
	}
	inline BinnedVector<T>& addCol(const BinnedVector<T>& col) { 
		if(col.bin() != m_bin)
			throw std::runtime_error("Incompatible bins.");
		m_y.push_back(col);
		return m_y.back();
	}

	std::ostream& printY(std::ostream& o, int index) const {
		for(auto a = m_y.begin(); a != m_y.end(); ++a)
		{
			if(a->size() > index)
				o << '\t' << (*a)()[index];
			else
				o << "\tNAN";
		}
		return o;
	}
	std::ostream& print(std::ostream& o, int index) const {
		o << index*m_bin;
		return printY(o, index);
	}
	std::ostream& print(std::ostream& o, int index, double xoffset) const {
		o << xoffset+index*m_bin;
		return printY(o, index);
	}
	std::ostream& print(std::ostream& o, double xoffset) const {
		const int tsize = size();
		for(int a = 0; a < tsize; ++a)
			print(o, a, xoffset) << '\n';
		return o;
	}

	int size() const {
		size_t n = 0;
		for(auto a = m_y.begin(); a != m_y.end(); ++a)
			n = std::max(n, a->size());
		return n;
	}
	inline void resize(double size) {
		for(auto a = m_y.begin(); a != m_y.end(); ++a)
			a->resize(size);
	}
	inline void resize(double size, const T& value) {
		for(auto a = m_y.begin(); a != m_y.end(); ++a)
			a->resize(size, value);
	}
	inline int cols() const { return m_y.size();}
	inline BinnedVector<T>& col(int n) {return m_y[n];}

};

class BinnedVectorWindow
{
	protected:
	double m_min, m_size;

	public:

		BinnedVectorWindow(double min, double size)
			: m_min(min), m_size(size)
		{}

	template<class T>
	inline T& operator() (BinnedVector<T>& vec, double i) const  {return vec[i-m_min];}
	template<class T>
	inline const T& operator() (const BinnedVector<T>& vec, double i) const {return vec[i-m_min];}

	template<class T, template<class Tt> class BinnedContainer>
	inline void resize(BinnedContainer<T>& vec) {
		vec.resize(m_size);
	}
	template<class T, template<class Tt> class BinnedContainer>
	inline void resize(BinnedContainer<T>& vec, const T& value) {
		vec.resize(m_size+vec.bin(), value);
	}

	inline double min() const {return m_min;}
	inline double size() const {return m_size;}
};

class BinnedVectorWindowCentered : public BinnedVectorWindow
{
	double m_halfBin;

	public:

		BinnedVectorWindowCentered(double min, double max, double bin)
			: BinnedVectorWindow(min-bin/2., max+bin/2.-min), m_halfBin(bin/2.)
		{}

		inline double max() const {return m_min+m_size;}
		inline double min() const {return m_min+m_halfBin;}
};

template<class T>
std::ostream& operator<< (std::ostream& o, const MultiColumnBinnedVector<T>& dat)
{
	const int size = dat.size();
	for(int a = 0; a < size; ++a)
		dat.print(o, a) << '\n';
	return o;
}

template<class T>
std::ostream& operator<< (std::ostream& o, const BinnedVector<T>& dat)
{
	double x = 0;
	for(int a = 0; a < dat.size(); ++a, x += dat.bin())
		o << x << '\t' << dat()[a] << '\n';
	return o;
}
#endif
