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

#ifndef KMC_COMBINE_TABLES_H
#define KMC_COMBINE_TABLES_H

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <list>
#include <vector>
#include <cmath>

#include "binnedVector.h"

namespace Kmc
{

class OStreamSwitch
{
	bool m_on;
	std::ostream& m_stream;
	std::string m_prefix;
	static std::ostream null;

	public:

		OStreamSwitch(const std::string& prefix = std::string(), std::ostream& stream = std::cerr, bool on = true)
			: m_on(on), m_stream(stream), m_prefix(prefix)
	{}

	void on() {m_on = true;}
	void off() {m_on = false;}

	std::ostream& operator() () {
		if(m_on)
			return m_stream << m_prefix;
		return null;
	}
};
extern OStreamSwitch WARN;

class TableFileView
{
	public:
	class Pair {
		double m_x, m_sig;
		std::vector<double> m_y;
		double m_n;
		std::string m_meta;

		public:

			Pair() {}
			Pair(double x, int ny = 0)
				: m_x(x), m_y(ny), m_n(0), m_sig(0.) {}
			Pair(double x, const std::vector<double>& y, double n = 0., double sig = 0)
				: m_x(x), m_y(y), m_n(n), m_sig(sig) {}
			Pair(const Pair& other)
				: m_x(other.m_x), m_sig(other.m_sig), m_y(other.m_y), m_n(other.m_n) {}

		double x() const {return m_x;}
		double& rx() {return m_x;}
		const std::vector<double>& y() const {return m_y;}
		double y(int i) const {return m_y.at(i);}
		std::vector<double>& ry() {return m_y;}
		double& ry(int i) {return m_y.at(i);}
		double n() const {return m_n;}
		double& rn() {return m_n;}
		const std::string& meta() const {return m_meta;}
		std::string& rmeta() {return m_meta;}
		double sig() const {return m_sig;}
		double& rsig() {return m_sig;}
		double var() const {return m_sig*m_sig;}
		double sumSq(int i) const {return y(i)*y(i)*n() + var()*(n()-1);}
		int ny() const {return m_y.size();}
	};

	private:
	int m_idxX, m_idxN, m_idxSigma, m_maxIdx, m_idxMeta;
	std::vector<int> m_idxY;

	std::string m_name, m_lastLine;
	std::ifstream m_file;

	int m_line;
	Pair m_currentPair;

	struct Line;

	public:

	TableFileView(int idxX, std::list<int> &idxY, const char* name, int idxN = -1, int m_idxSigma = -1, int idxMeta = -1);
	TableFileView(int idxX, int idxY, const char* name, int idxN = -1, int m_idxSigma = -1, int idxMeta = -1);
	/*! No y column, consequently no N and sigma. */
	TableFileView(int idxX, const char* name, int idxMeta = -1);

	bool good() const {return m_file.good();}
	double x() const {return m_currentPair.x();}
	int idxX() const {return m_idxX;}
	const std::vector<double>& y() const {return m_currentPair.y();}
	double y(int i) const {return m_currentPair.y(i);}
	int idxY(int i) const {return m_idxY.at(i);}
	int nY() const {return m_currentPair.ny();}
	double n() const {return m_currentPair.n();}
	int idxN() const {return m_idxN;}
	double sig() const {return m_currentPair.sig();}
	int idxSigma() const {return m_idxSigma;}
	const std::string& meta() const {return m_currentPair.meta();}
	int idxMeta() const {return m_idxMeta;}
	int line() const {return m_line;}
	const Pair& currentPair() const {return m_currentPair;}
	const std::string& name () const {return m_name;}
	const std::string& lastLine () const {return m_lastLine;}
	bool open();
	void close() {m_file.close();}
	TableFileView& nextLine();

	/*! Read complete contents of associated file.
	 *
	 * Will close and reopen file already opened
	 * \param xarr vector for x
	 * \param yarr vector of vectors for y cols
	 * \return success
	 */
	bool read(std::vector<double>& xarr, std::vector<std::vector<double> >& yarr);
};

bool parseIdxYList(std::list<int>& idxY, const std::string& list);

class Sample
{
	protected:
	long long m_x;
	double m_sum, m_sumSq, m_n;
	const TableFileView* m_file;

	public:

		Sample(long long x, const TableFileView* file = 0)
			: m_x(x), m_n(0.), m_sum(0.), m_sumSq(0.), m_file(file) {}

	long long x() const {return m_x;}
	double  n() const {return m_n;}
	double sum() const {return m_sum;}
	double sumSq() const {return m_sumSq;}
	double avg() const {return m_sum/m_n;}
	double avgSq() const {return m_sumSq/m_n;}
	double stdev() const {return sqrt((m_sumSq - m_sum*m_sum/m_n)/(m_n-1));}
	const TableFileView* file() const {return m_file;}

	void setFile(const TableFileView* file) {m_file = file;}

	Sample& operator+=(double y);
	Sample& addWeighted(double y, double weight);
	Sample& operator+=(const TableFileView& view);
	Sample& addWeighted(const TableFileView& view, double weight);
};

class CacheSample : public Sample
{
	std::vector<double> m_values;
	double m_max, m_min;

	public:

		CacheSample(long long x, const TableFileView* file = 0, unsigned int initialNValues = 2000)
			: Sample(x, file), m_min(INFINITY),  m_max(-INFINITY) {
			m_values.reserve(initialNValues);
		}

	double min() const {return m_min;}
	double max() const {return m_max;}
	const std::vector<double>& values() const {return m_values;}
	bool histGood() const {return n()>0 && m_min < m_max;}
	unsigned int nbinsFromWidth(double width) const {return ceil((m_max-m_min)/width);}

	std::pair<BinnedVector<unsigned int>, BinnedVectorWindow> makeHist(unsigned int nbins);
	void printHist(std::ostream& o, unsigned int nbins);

	CacheSample& operator+=(double y);
	CacheSample& addWeighted(double y, double weight);
	CacheSample& operator+=(const TableFileView& view);
	CacheSample& addWeighted(const TableFileView& view, double weight);
};

}

std::ostream& operator<<(std::ostream& o, const Kmc::TableFileView& a);
std::ostream& operator<<(std::ostream& o, const Kmc::Sample& s);
std::ostream& operator<<(std::ostream& o, const Kmc::CacheSample& s);

#endif
