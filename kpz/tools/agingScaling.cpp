#include <combineTables.h>
#include <arghelper.h>

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

using namespace Kmc;

inline void eraseFromVec(std::vector<double>& vec, unsigned int min, unsigned int sup)
{
	vec.erase(vec.begin()+min, vec.begin()+sup);
}

int main(int argc, char* argv[])
{
	if(argc < 5)
	{
		std::cerr << "Usage aging (ref) (s1) [-n norm] (interpol) (s2) [-n norm] [options]\n";
		return 1;
	}

	std::vector<double> vecRefX;
	std::vector<std::vector<double> > vecY;
	double s[2];
	unsigned int interpolNameIdx = 3;
	double xMin = 1.00001, xMax = INFINITY;
	double eMin = -0.35;
	double eStep = 0.01;
	bool plInt = false;
	unsigned int eStepN = 15;

	if(!getArg(s[0], 2, argc, argv))
		return 1;

	{
		TableFileView tfvRef(0, 1, argv[1]);
		tfvRef.read(vecRefX, vecY);
		for(auto& a : vecRefX)
			a /= s[0];
		if(!strcmp(argv[3], "-n"))
		{
			double tmp;
			if(!getArg(tmp, 4, argc, argv))
			{
				std::cerr << "[EE] expected double value for -n";
				return 1;
			}
			for(auto& a : vecY[0])
				a /= tmp;
			interpolNameIdx += 2;
		}
	}
	if(!getArg(s[1], interpolNameIdx+1, argc, argv))
		return 1;
	{
		std::vector<double> vecInterpolX;
		std::vector<std::vector<double> > vecInterpolY;
		TableFileView tfvRef(0, 1, argv[interpolNameIdx]);
		tfvRef.read(vecInterpolX, vecInterpolY);
		vecY.push_back(vecY[0]);
		for(auto& a : vecInterpolX)
			a /= s[1];

		if(interpolNameIdx+2 < argc && !strcmp(argv[interpolNameIdx+2], "-n"))
		{
			double tmp;
			if(!getArg(tmp, interpolNameIdx+3, argc, argv))
			{
				std::cerr << "[EE] expected double value for -n";
				return 1;
			}
			for(auto& a : vecInterpolY[0])
				a /= tmp;
			interpolNameIdx += 2;
		}

		for(int i = interpolNameIdx+2; i < argc; ++i)
		{
			if(!strcmp(argv[i], "--xMin"))
			{
				if(!getArg(xMin, ++i, argc, argv))
					return 1;
			}
			else if(!strcmp(argv[i], "--xMax"))
			{
				if(!getArg(xMax, ++i, argc, argv))
					return 1;
			}
			else if(!strcmp(argv[i], "--eMin"))
			{
				if(!getArg(eMin, ++i, argc, argv))
					return 1;
			}
			else if(!strcmp(argv[i], "--plInt"))
			{
				plInt = true;
			}
			else if(!strcmp(argv[i], "--eIter"))
			{
				if(!getArg(eStep, ++i, argc, argv))
					return 1;
				if(!getArg(eStepN, ++i, argc, argv))
					return 1;
			}
			else
			{
				std::cerr <<"[EE] Unknown option: " << argv[i] << '\n';
				return 1;
			}
		}

		unsigned int refSup = vecRefX.size(), refMin = 0;
		unsigned int interMin = 0;
		for(; interMin < vecInterpolX.size(); ++interMin)
			if(vecInterpolX[interMin] >= xMin)
				break;
		for(int inter = 0; inter < vecInterpolX.size(); ++inter)
		{
			if(vecRefX[inter] <= 1)
				continue;
			if(vecRefX[inter] >= vecInterpolX[interMin])
			{
				refMin = inter;
				break;
			}
		}

		std::ofstream out("agingScaling_inter.dat");

		for(int ref = refMin, inter = interMin; ref < vecRefX.size(); ++ref)
		{
			if(vecRefX[ref] > xMax)
			{
				refSup = ref;
				break;
			}
			while(inter < vecInterpolX.size()-1 && vecInterpolX[inter+1] < vecRefX[ref])
				++inter;
			if(inter >= vecInterpolX.size())
			{
				refSup = ref;
				break;
			}
			if(plInt)
				vecY[1][ref] = exp( log(vecInterpolY[0][inter]) + (log(vecInterpolY[0][inter+1])-log(vecInterpolY[0][inter]))
					* (log(vecRefX[ref])-log(vecInterpolX[inter])) / (log(vecInterpolX[inter+1])-log(vecInterpolX[inter])));
			else
				vecY[1][ref] = vecInterpolY[0][inter] + (vecInterpolY[0][inter+1]-vecInterpolY[0][inter])
					* (vecRefX[ref]-vecInterpolX[inter]) / (vecInterpolX[inter+1]-vecInterpolX[inter]);
			if(vecY[0][ref] <= 0. || vecY[1][ref] <= 0. || isnan(vecY[1][ref]))
			{
				refSup = ref;
				break;
			}
			out << vecRefX[ref] << '\t' << vecY[0][ref] << '\t' << vecY[1][ref] << '\n';
		}
		if(refMin > 0)
		{
			eraseFromVec(vecRefX, 0u, refMin);
			eraseFromVec(vecY[0], 0u, refMin);
			eraseFromVec(vecY[1], 0u, refMin);
			refSup -= refMin;
		}
		if(refSup > 0)
		{
			eraseFromVec(vecRefX, refSup, vecRefX.size());
			for(int a = 0; a < 2; ++a)
				eraseFromVec(vecY[a], refSup, vecY[a].size());
		}
	}

	double e = eMin;
	const double sRatio = log(s[0]/s[1]);
	{
		double sum = 0, sum2 = 0, weights = 0;
		for(int a = 0; a < vecRefX.size(); ++a)
		{
			const double d = log(vecY[0][a] / vecY[1][a])/sRatio;
			const double w = vecY[0][a] * vecY[1][a];
			sum += d*w;
			sum2 += d*d*w;
			weights += w;
		}
		sum /= weights;
		sum2 /= weights;
		double stddev = sqrt(sum2 - sum*sum);
		std::cout << "#e=" << sum << '\t' << stddev << '\n';// << stddev/sqrt(vecRefX.size()) << '\n';
	}
	for(int a = 0; a < eStepN; ++a, e+=eStep)
	{
		std::cout << e << '\t';
		
		double chi2 = 0, weights = 0;
		for(int a = 0; a < vecRefX.size(); ++a)
		{
			const double d = log(vecY[0][a] / vecY[1][a])/sRatio - e;
			const double w = vecY[0][a] * vecY[1][a];
			chi2 += d*d*w;
			weights += w;
		}
		std::cout << chi2/weights << '\n';
	}

	return 0;
}
