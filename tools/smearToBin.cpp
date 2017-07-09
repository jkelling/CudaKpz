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

#include <combineTables.h>
#include <binnedVector.h>
#include <binnedVector.h>
#include <arghelper.h>

#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <list>

using namespace Kmc;

class Gaussian
{
	double m_a, m_b, m_cSq;

	public:

		Gaussian() {}
		Gaussian(double a, double b, double c)
			{set(a,b,c);}
	/*! Construct normalized gaussian. */
		Gaussian(double mu, double sig)
			{setNormalized(mu, sig);}

	void set(double a, double b, double c) {
		m_a = a;
		m_b = b;
		m_cSq = c*c;
	}
	void setNormalized(double mu, double sig) {
		set(1./(sig*sqrt(2*M_PI)), mu, sig);
	}

	double operator() (double x) const {
		const double numeratorRoot = x-m_b;
		return m_a*exp(-numeratorRoot*numeratorRoot/2/m_cSq);
	}

	double mu() const {return m_b;}
	double sig() const {return sqrt(m_cSq);}
};

template<class Smearing>
void smearToBin(const Smearing& smearing, BinnedVector<double>& outData, const BinnedVectorWindowCentered& window, const std::vector<double>& x, const std::vector<double>& y)
{
	const double max = window.max();
	for(double a = window.min(); a <= max; a += outData.bin())
	{
		double sum = (x[1] - x[0]) * y[0] * smearing(x[0]);
		for(int i = 1; i < x.size()-1; ++i)
		{
			sum += (x[i+1] - x[i-1])/2. * y[i] * smearing(a-x[i]);
		}
		sum += (x[x.size()-1] - x[x.size()-1]) * y[x.size()-1] * smearing(x[x.size()-1]);
		window(outData, a) += sum;
	}
}

int main(int argc, char* argv[])
{
	double xOutMin = 0, xOutMax = 10, xOutBin = .1;
	int idxX = 0;
	std::list<int> idxY{1};
	int fnameIdx = 0;
	Gaussian gauss(0.,0.01);

	for(int i = 1; i < argc; ++i)
	{
		if(!strcmp("-x", argv[i]))
		{
			if(!getArg(idxX, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("-y", argv[i]))
		{
			std::string tmp;
			if(!getArg(tmp, ++i, argc, argv))
				return 1;
			if(!parseIdxYList(idxY, tmp))
				return 1;
		}
		else if(!strcmp("-o", argv[i]))
		{
			if(!getArg(xOutMin, ++i, argc, argv))
				return 1;
			if(!getArg(xOutMax, ++i, argc, argv))
				return 1;
			if(!getArg(xOutBin, ++i, argc, argv))
				return 1;
		}
		else if(!strcmp("-gauss", argv[i]))
		{
			double mu = 0.,  sig;
			if(!getArg(sig, ++i, argc, argv))
				return 1;
			if(!getArg(mu, ++i, argc, argv))
				--i;
			gauss.setNormalized(mu, sig);
		}
		else
		{
			if(fnameIdx == 0)
				fnameIdx = i;
			else
			{
				std::cerr << "[EE] Unknown argument: " << argv[i] << '\n';
				return 1;
			}
		}
	}
	if(fnameIdx == 0)
	{
		std::cerr << "[EE] No indpufile given\n";
		return 1;
	}

	std::vector<double> x;
	std::vector<std::vector<double> > y;
	{
		TableFileView file(idxX, idxY, argv[fnameIdx]);
		if(!file.read(x, y))
		{
			std::cerr << "[EE] failed to read file: " << argv[fnameIdx] << '\n';
			return 1;
		}
	}
	if(x.size()<3)
	{
		std::cerr << "[EE] Too few points in input (<3): " << x.size() << '\n';
		return 1;
	}
	
	MultiColumnBinnedVector<double> outData(xOutBin, idxY.size());
	BinnedVectorWindowCentered window(xOutMin, xOutMax, outData.bin());
	window.resize(outData, 0.);

	for(int a = 0; a < y.size(); ++a)
	{
		smearToBin(gauss, outData.col(a), window, x, y[a]);
	}

	std::cout << "#smearing: gauss ,sig= " << gauss.sig() << " ,mu= " << gauss.mu() << '\n';
	outData.print(std::cout, window.min());

	return 0;
}
