/***************************************************************************
*   Copyright 2011 - 2016 Jeffrey Kelling <j.kelling@hzdr.de>
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

#include "combineTables.h"

#include <cmath>
#include <algorithm>

std::ostream Kmc::OStreamSwitch::null(0);
Kmc::OStreamSwitch Kmc::WARN("[WW] ", std::cerr, false);

struct Kmc::TableFileView::Line
{
	int col = 0, yCol = 0;
	std::istringstream in;
	TableFileView& view;

	bool skip() {
		if(!std::getline(view.m_file, view.m_lastLine))
			return false;
		in.clear();
		in.str("");
		in.str(view.m_lastLine);
		++view.m_line;
		col = yCol = 0;
		return true;
	}

		Line(TableFileView& view) : view(view) {
		skip();
	}
};

Kmc::TableFileView::TableFileView(int idxX, std::list<int> &idxY, const char* name, int idxN, int idxSigma, int idxMeta)
	: m_idxX(idxX), m_name(name), m_line(0), m_idxN(idxN), m_idxSigma(idxSigma), m_idxMeta(idxMeta)
{
	if(!m_idxY.empty())
	{
		m_maxIdx = std::max(m_idxX, m_idxMeta);
		m_idxN = m_idxSigma = -1;
		return;
	}

	idxY.sort();
	idxY.unique();
	auto a = idxY.begin();
	for(; a != idxY.end() && *a < m_idxX; ++a);
	if(a != idxY.end() && *a == m_idxX)
		idxY.erase(a);

	const int nY = idxY.size();
	m_idxY.resize(nY);
	m_currentPair.ry().resize(nY);
	int i = 0;
	for(a = idxY.begin(); a != idxY.end(); ++i, ++a)
	{
		m_idxY[i] = *a;
	}

	// assume single samples if not n-col is available
	if(m_idxN < 0)
		m_currentPair.rn() = 1;

	m_maxIdx = std::max(std::max(m_idxX, m_idxN), std::max(m_idxMeta, *std::max_element(m_idxY.begin(), m_idxY.end())));
}

Kmc::TableFileView::TableFileView(int idxX, int idxY, const char* name, int idxN, int idxSigma, int idxMeta)
	: m_idxX(idxX), m_name(name), m_line(0), m_idxN(idxN), m_idxSigma(idxSigma), m_idxMeta(idxMeta)
{
	const int nY = 1;
	m_idxY.resize(nY);
	m_idxY[0] = idxY;
	m_currentPair.ry().resize(nY);
	int i = 0;

	// assume single samples if not n-col is available
	if(m_idxN < 0)
		m_currentPair.rn() = 1;

	m_maxIdx = std::max(std::max(m_idxX, m_idxN), std::max(m_idxMeta, m_idxY[0]));
}

Kmc::TableFileView::TableFileView(int idxX, const char* name, int idxMeta)
	: m_idxX(idxX), m_name(name), m_line(0), m_idxN(-1), m_idxSigma(-1), m_idxMeta(idxMeta)
{
	m_maxIdx = std::max(m_idxX, m_idxMeta);
}

bool Kmc::TableFileView::open()
{
	m_file.open(m_name);
	nextLine();
	return good();
}

Kmc::TableFileView& Kmc::TableFileView::nextLine()
{
	char c;
	std::string tmp;
	while(good())
	{
		Line line(*this);
	while(line.in >> c)
	{
		if(c == '#')
		{
			if(!line.skip())
				return *this;
			continue;
		}
		line.in.putback(c);
		if(line.col == m_idxX)
		{
			if(!(line.in >> m_currentPair.rx()))
			{
				WARN() << m_name << ": skipping defective line " << m_line << '\n';
				if(!line.skip())
					return *this;
				continue;
			}
		}
		else if(line.col == m_idxY[line.yCol])
		{
			if(!(line.in >> m_currentPair.ry(line.yCol)))
			{
				WARN() << m_name << ": skipping defective line " << m_line << '\n';
				if(!line.skip())
					return *this;
				continue;
			}
			++line.yCol;
		}
		else if(m_idxN > 0 && line.col == m_idxN)
		{
			if(!(line.in >> m_currentPair.rn()))
			{
				WARN() << m_name << ": skipping defective line " << m_line << '\n';
				if(!line.skip())
					return *this;
				continue;
			}
		}
		else if(m_idxSigma > 0 && line.col == m_idxSigma)
		{
			if(!(line.in >> m_currentPair.rsig()))
			{
				WARN() << m_name << ": skipping defective line " << m_line << '\n';
				if(!line.skip())
					return *this;
				continue;
			}
		}
		else if(m_idxMeta > 0 && line.col == m_idxMeta)
		{ //! \todo parse ""
			if(!(line.in >> m_currentPair.rmeta()))
			{
				WARN() << m_name << ": skipping defective line " << m_line << '\n';
				if(!line.skip())
					return *this;
				continue;
			}
		}
		else
		{
			std::string tmp;
			if(!(line.in >> tmp))
			{
				WARN() << m_name << ": skipping defective line " << m_line << '\n';
				if(!line.skip())
					return *this;
			}
		}
		++line.col;

		if(line.col > m_maxIdx)
			return *this;
	}
		if(line.col > m_maxIdx)
			break;
	}
	return  *this;
}

bool Kmc::TableFileView::read(std::vector<double>& xarr, std::vector<std::vector<double> >& yarr)
{
	if(m_file.is_open())
		m_file.close();
	if(!open())
		return false;

	yarr.resize(nY());
	xarr.clear();
	int capacity = 100;
	if(xarr.capacity() >= capacity)
		capacity = xarr.capacity();
	else
		xarr.reserve(capacity);
	for(auto& a : yarr)
		a.reserve(capacity);

	while(good())
	{
		xarr.push_back(x());
		for(int a = 0; a < nY(); ++a)
			yarr[a].push_back(y(a));
		if(xarr.size() == capacity)
		{
			capacity <<=1;
			xarr.reserve(capacity);
			for(auto& a : yarr)
				a.reserve(capacity);
		}
		nextLine();
	}
	xarr.shrink_to_fit();
	for(auto& a : yarr)
		a.shrink_to_fit();
	return true;
}

bool Kmc::parseIdxYList(std::list<int>& idxY, const std::string& list)
{
	std::istringstream in(list);
	idxY.clear();
	unsigned int tmp;
	while(!(in.eof() || in.bad()))
	{
		char c = 0;
		in.clear();
		in >> c;
		switch(c)
		{
			case ',':
				break;
			case '-':
				if((in >> tmp).fail())
				{
					std::cerr << "[EE] Expected integer after '-'\n";
					return false;
				}
				for(int a = idxY.back() + 1; a <= tmp; ++a)
				{
					idxY.push_back(a);
				}
				break;
			default:
				if(isdigit(c))
					in.putback(c);
				else
				{
					std::cerr << "[EE] Invalid charakter in idx list: " << c << '\n';
					return false;
				}
		}
		if (!(in >> tmp).fail())
			idxY.push_back(tmp);
	}
	return true;
}

Kmc::Sample& Kmc::Sample::operator+=(double y)
{
	++m_n;
	m_sum += y;
	m_sumSq += y*y;
	return *this;
}

Kmc::Sample& Kmc::Sample::addWeighted(double y, double weight)
{
	const double sumSq_v = y*y*weight;
	const double n_new = m_n+weight;
	m_sum += y*weight;
	const double sqSum = m_sum*m_sum;
	m_sumSq = ((m_sumSq + sumSq_v)*n_new - sqSum)*(n_new-1)
			/ (n_new - m_n - weight)
			+ sqSum/n_new;
	m_n = n_new;
	return *this;
}

Kmc::Sample& Kmc::Sample::operator+=(const TableFileView& view)
{
	static const int I = 0;
	if(view.n() == 1)
	{
		for(int y = 0; y < view.nY(); ++y)
		{
			(*this) += view.y(y);
		}
		return *this;
	}
	m_n += view.n();
	m_sum += view.y(I)*view.n();
	m_sumSq += view.currentPair().sumSq(I);
	return *this;
}

Kmc::Sample& Kmc::Sample::addWeighted(const TableFileView& view, double weight)
{
	static const int I = 0;
	if(m_n > 0)
	{
		const double n_v = view.n()*weight;
		const double sumSq_v = view.currentPair().sumSq(I)*weight;
		const double n_new = m_n+n_v;
		m_sum += view.y(I)*n_v;
		const double sqSum = m_sum*m_sum;
		m_sumSq = ((m_sumSq + sumSq_v)*n_new - sqSum)*(n_new-1)
				/ (n_new*n_new - m_n - n_v*weight)
				+ sqSum/n_new;
		m_n = n_new;
	}
	else
	{
		m_n = view.n()*weight;
		m_sumSq = view.currentPair().sumSq(I)*weight;
		m_sum += view.y(I)*m_n;
	}
	return *this;
}

Kmc::CacheSample& Kmc::CacheSample::operator+=(double y)
{
	Sample::operator+=(y);
	if(m_min > y)
		m_min = y;
	else if(m_max < y)
		m_max = y;
	m_values.push_back(y);
	if(m_values.capacity() <= m_values.size())
		m_values.reserve(m_values.size()+2000);
	return *this;
}

Kmc::CacheSample& Kmc::CacheSample::addWeighted(double y, double weight)
{
	Sample::addWeighted(y,weight);
	(*this) += y;
	return *this;
}

Kmc::CacheSample& Kmc::CacheSample::operator+=(const TableFileView& view)
{
	static const int I = 0;
	if(view.n() == 1)
	{
		for(int y = 0; y < view.nY(); ++y)
		{
			(*this) += view.y(y);
		}
		return *this;
	}
	std::cerr << "[WW][CacheSample::+=] Will not add multiple copies.\n";
	m_n += view.n();
	m_sum += view.y(I)*view.n();
	m_sumSq += view.currentPair().sumSq(I);
	return (*this) += view.y(I);
}

Kmc::CacheSample& Kmc::CacheSample::addWeighted(const TableFileView& view, double weight)
{
	std::cerr << "[WW][CacheSample::+=] Will not add weighted sample.\n";
	return (*this)+=view;
}

std::pair<BinnedVector<unsigned int>, BinnedVectorWindow> Kmc::CacheSample::makeHist(unsigned int nbins)
{
	BinnedVector<unsigned int> hist((m_max-m_min)/(double)nbins);
	BinnedVectorWindow histWindow {m_min, m_max-m_min};
	if(!histGood())
	{
		std::cerr << "[EE][CacheSample::makeHist] Insufficient data.\n";
		return std::move(std::make_pair(std::move(hist), std::move(histWindow)));
	}
	histWindow.resize(hist, 0u);

	for(auto a : m_values)
		++histWindow(hist, a);

	return std::move(std::make_pair(std::move(hist), std::move(histWindow)));
}

void Kmc::CacheSample::printHist(std::ostream& o, unsigned int nbins)
{
	auto h = makeHist(nbins);

	unsigned int count = 0;
	for(int a = 0; a < h.first.size(); ++a)
		count += h.first.at(a);

	o << "# N=" << count << "\n# avg=" << avg() << '\n';
	for(int a = 0; a < h.first.size(); ++a)
	{
		const double val = h.first.at(a)/(double)count/h.first.bin();
		o << h.second.min()+a*h.first.bin() << '\t' << val << '\t';
		o << (h.second.min()+a*h.first.bin())/avg() << '\t' << val*avg() << '\n';
	}
}

std::ostream& operator<<(std::ostream& o, const Kmc::TableFileView& a)
{
	o << "TableFileView " << a.name () << "\n\tx = " << a.idxX() << "\n\ty =";
	for(int i = 0; i < a.nY(); ++i)
		o << ' ' << a.idxY(i);
	return o << '\n'; 
}

std::ostream& operator<<(std::ostream& o, const Kmc::Sample& s)
{
	o << s.x() << '\t';
	if(s.n() > 1)
	{
		const double stdev = s.stdev();
		o << s.avg() << '\t' << stdev << '\t'
			<< stdev/sqrt((double)s.n()) << '\t' << s.n();
	}
	else
	{
		o << s.sum() << "\t0\t0\t" << s.n() << "\t#";
		if(s.file())
			o << ' ' << s.file()->name();
	}
	return o;
}

std::ostream& operator<<(std::ostream& o, const Kmc::CacheSample& s)
{
	o << (const Kmc::Sample&)s;
	if(s.n() > 1)
	{
		double k[4] = {0.,0.,0.,0.};
		for(auto v : s.values())
		{
			const double d = v-s.avg();
			double acc = d;
			for(int i = 0; i < 4; ++i, acc *= d)
				k[i] += acc;
		}
		for(int i = 0; i < 4; ++i)
			k[i] /= s.n();

		o << '\t' << (k[2] / pow(k[1],1.5)) << '\t' << (k[3] / (k[1] * k[1]) - 3.0);
	}
	return o;
}
