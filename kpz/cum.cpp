/***************************************************************************
*   Copyright 2011 Geza Odor <odor@mfa.kfki.hu>
*                  Research Centre for Natural Sciences
*                  Hungarian Academy of Sciences
*   Copyright 2011 Jeffrey Kelling <j.kelling@hzdr.de>
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

#include <arghelper.h>

#include <cstdio>
#include <cmath>
#include <cstring>

#include <vector>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <sstream>

#define tmax 1000000
#define Hmax 200

int j,k,ca,navar,np,ndum,h1i;
long double w1,w2,w3,w4,k1,k2,k3,k4,hmean,htot,cin;
float p,z;
double o1,o2,o3,o4;
long double natl;
FILE *fopen(),*fp1, *fp0;
char *fst;

std::istream& get(std::istream& i, int& tin, long double& cin, const int col)
{
  while(!(i >> tin))
  {
    if(i.rdstate() == std::ios::failbit)
    {
      i.clear();
      i.ignore(100000, '\n');
    }
    else
      return i;
  }
  std::string tmp;
  std::getline(i, tmp);
  std::istringstream is(tmp);
  for(int a = 0; a < col; ++a)
  { 
    if(!(is >> cin))
    {
      i.setstate(std::ios::failbit);
      break;
    }
  }
  return i;
}

main (int argc, char* argv[]) {

  bool lazy = false;
  int col = 1, nbin = 0;
  long double hmin = 0, hmax = 0;
  std::list<std::string> filelist;
  for(int i = 1; i < argc; ++i)
  {
    if(!strcmp(argv[i],"--lazy"))
      lazy = true;
    else if(!strcmp(argv[i],"-h"))
    {
      if(!getArg(nbin, ++i, argc, argv))
        return 1;
      if(!getArg(hmin, ++i, argc, argv))
        return 1;
      if(!getArg(hmax, ++i, argc, argv))
        return 1;
    }
    else if(!strcmp(argv[i], "--filelist"))
    {
      std::string tmp;
      if(!getArg(tmp, ++i, argc, argv))
        return 1;
      std::ifstream i(tmp.c_str());
      if(!i)
      {
        std::cerr << "could not obtain filelist from " << tmp << '\n';
        return 1;
      }
      while(std::getline(i, tmp))
        filelist.push_back(tmp);
    }
    else if(!strcmp(argv[i], "-f"))
    {
      std::string tmp;
      if(!getArg(tmp, ++i, argc, argv))
        return 1;
      filelist.push_back(tmp);
    }
    else if(!strcmp(argv[i], "-a") || !strcmp(argv[i], "--avertime"))
    {
      if(!getArg(navar, ++i, argc, argv))
        return 1;
    }
    else if(!strcmp(argv[i], "-c"))
    {
      if(!getArg(col, ++i, argc, argv))
        return 1;
    }
    else if(!strcmp(argv[i], "--help"))
    {
      std::cout << "usage: cum [options] [datafiles]\n"
        "options:\n"
        " -a / --avertime time\n"
        " -h nbin hmin hmax\n"
        " -c col\n"
        " -f datafile\n"
        " --filelist filelistfile\n"
        " --lazy                     don't demand times to match across files\n"
        " --help                     show this help and exit\n";
      return 0;
    }
    else
    {
      filelist.push_back(argv[i]);
    }
  }
  
  std::vector<long double> W1, W2, W3, W4, t, hist;
  if(lazy)
  {
	  W1.resize(tmax);
	  W2.resize(tmax);
	  W3.resize(tmax);
	  W4.resize(tmax);
	  t.resize(tmax);
  }
  else
  {
	  W1.reserve(tmax);
	  W2.reserve(tmax);
	  W3.reserve(tmax);
	  W4.reserve(tmax);
	  t.reserve(tmax);
  }
  htot = 0; natl=0; hmean=0;
  hist.resize(nbin, 0.); // allocate memory for hist and init with 0.
  const long double dh = (hmax-hmin)/nbin;

  int fnum = 0;
  for (; filelist.size(); filelist.pop_front()) {
    std::ifstream i (filelist.front().c_str());
    if(!i)
    {
      std::cerr << "Error: could not open " << filelist.front() << '\n';
      continue;
    }
    int tin;
    if(!fnum && !lazy)
    {
      while(get(i, tin, cin, col))
      {
        w1 = cin - 0.25;
        t.push_back(tin);
        W1.push_back(w1);
        W2.push_back(w1*w1);
        W3.push_back(w1*w1*w1);
        W4.push_back(w1*w1*w1*w1);

        if(tin >= navar) {
          const int h = (w1 - hmin)/dh;
          if(h >= 0 && h < hist.size())
            ++hist[h];
        }
     }
     ++fnum;
   }
   else
   {
     int a, skip = 0;
     for(a = 0; a < t.size() && get(i, tin, cin, col); )
     {
        w1 = cin - 0.25;
        if(!lazy && t[a] != tin)
        {
          std::cerr << "Fatal: input file " << filelist.front() << " contains incompatible data.\n";
          std::cerr << a << ' ' << t[a] << ' ' << tin << '\n';
          return 1;
        }
        else if(lazy)
        {
          if(tin < t[a])
          {
            ++skip;
            continue;
          }
        }
        W1[a] += w1;
        W2[a] += w1*w1;
        W3[a] += w1*w1*w1;
        W4[a] += w1*w1*w1*w1;

        if(tin >= navar) {
          const int h = (w1 - hmin)/dh;
          if(h >= 0 && h < hist.size())
            ++hist[h];
        }
        ++a;
    }
    if(skip > 0)
    {
      std::cerr << "Warning: " << skip << " values in input file " << filelist.front() << " were skipped\n";
    }
	if(!lazy)
	{
    if(a < t.size())
    {
      std::cerr << "Fatal: input file " << filelist.front() << " contained insufficient data.\n";
      return 1;
    }
    else if(get(i, tin, cin, col))
    {
      std::cerr << "Warning: input file " << filelist.front() << " was truncated at t = " << tin << '\n';
    }
	}
    ++fnum;
   }
  }
  int iavar = -1;
  for(int a = 0; a < t.size(); ++a)
  {
    W1[a] /= fnum;
    W2[a] /= fnum;
    W3[a] /= fnum;
    W4[a] /= fnum;
    if(iavar < 0 && t[a] >= navar)
      iavar = a;
  }
  if(iavar < 0)
    iavar = t.size();
  std::cout << "\nread " << fnum << " file(s) containing " << t.size() << " samples ( "
    << t.size()*fnum << " total ) " << (t.size()-iavar)*fnum << " relevant\n";

  fp1=fopen("aver.dat","w");
  for (int i = 1; i <= t.size(); i++) { 
    fprintf(fp1,"%Lf %f\n",t[i],sqrt(W2[i]));
  }
  fprintf(fp1,"stat = %d\n",fnum);
  fclose(fp1);

  w1 = 0; w2 = 0; w3 = 0; w4 = 0;
  for (int i = iavar; i <= t.size(); i++) {  
    w1 += W1[i];
    w2 += W2[i];
    w3 += W3[i];
    w4 += W4[i];
    natl++;
  }
  w1 /= natl; 
  w2 /= natl; 
  w3 /= natl; 
  w4 /= natl;

  for (int h = 0; h < nbin;  h++) htot += hist[h];
  htot = htot*dh;                    /* normalization of P(W^2) */
  for (int h = 0; h < nbin;  h++) hist[h] /= htot;

  for (int h = 0; h < nbin;  h++) hmean += (hist[h]*(hmin+h*dh)*dh);
  p = hmean;
  printf("<W^2> : %f\n",p); 

  k1 = w1;
  k2 = w2 - w1*w1;
  k3 = w3 - 3*w1*w2 + 2*w1*w1*w1;
  k4 = w4 - 4*w3*w1 - 3*w2*w2 + 12*w1*w1*w2 - 6*w1*w1*w1*w1;
  o1 = k1; o2 = k2; o3=k3; o4=k4;
  printf("%f %f %f %f %f\n",sqrt(fabs(k1)),o1,o2,o3,o4); 

  fp1=fopen("hist.dat","w");
  for (int h = 0; h < nbin;  h++) {
     if(hist[h]>0) {
       o1 = hmin+h*dh; o2 = hist[h];
      fprintf(fp1,"%e %e\n",o1,o2);
    }
  }
  fclose(fp1);

  fp1=fopen("psi.dat","w");
  for (int h = 0; h < nbin;  h++) {
    if(hist[h]>0) {
      o1 = (hmin+h*dh)/hmean; o2 = hist[h]*hmean;
      fprintf(fp1,"%e %e\n",o1,o2);
    } 
  }
  fclose(fp1);
 
  return 0;
}
