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

#ifndef KMC_TEMPLATE_META_HELPER_H
#define KMC_TEMPLATE_META_HELPER_H

// Darkest magick happening here...
/* Rekursive power, reduce depth trough bisection.*/
template<unsigned long long base, unsigned int exp, bool divByTwo> struct PowHelper;
template<unsigned long long base, unsigned int exp>
struct PowHelper<base, exp, true> {
	static const unsigned long long val = PowHelper<base, exp/2, exp/2%2==0>::val*PowHelper<base, exp/2, exp/2%2==0>::val;
};
template<unsigned long long base, unsigned int exp>
struct PowHelper<base, exp, false> {
	static const unsigned long long val = PowHelper<base, exp-1, (exp-1)%2==0>::val*base;
};
template<unsigned long long base, bool dummy>
struct PowHelper<base, 0, dummy> {
	static const unsigned long long val = 1;
};
template<unsigned long long base>
struct PowHelper<base, 1, false> {
	static const unsigned long long val = base;
};

// ... gets even darker here.
/* Rekursive sum, reduce depth trough bisection.
 * lowerBound excluded, upperBound included */
template< template<int>class step, int lowerBound, int upperBound, bool rangeDivByTwo, bool rangeNonEmpty> struct SumHelper;
template< template<int>class step, int lowerBound, int upperBound>
struct SumHelper<step, lowerBound, upperBound, true, true> {
	static const unsigned long long center = (lowerBound+upperBound)/2;
	static const bool next = ((center-lowerBound)%2==0);
	static const unsigned long long val = SumHelper<step, lowerBound, center, next, center!=lowerBound>::val + SumHelper<step, center, upperBound, next, center!=upperBound>::val;
};
template< template<int>class step, int lowerBound, int upperBound>
struct SumHelper<step, lowerBound, upperBound, false, true> {
	static const unsigned long long val = SumHelper<step, lowerBound, upperBound-1, (upperBound-1-lowerBound)%2==0, lowerBound!=(upperBound-1)>::val + step<upperBound>::val;
};
template< template<int>class step, int lowerBound, int upperBound>
struct SumHelper<step, lowerBound, upperBound, true, false> {
	static const unsigned long long val = 0;
};

template<unsigned int n>
struct Log2 {
	static const int val = 1 + Log2< (n>>1) >::val;
};
template<>
struct Log2<0>; // log(0) is not allowed
template<>
struct Log2<1> {
	static const int val = 0;
};
template<unsigned int n>
struct Log2RoundUp {
	static const int log = Log2<n>::val;
	static const int val = log + (int)((n<<(32-log))!=0);
};

template<bool cond, class If, class Else>
class TemplateCond;
template<class If, class Else>
struct TemplateCond<true, If, Else>{
	typedef If Result;
};
template<class If, class Else>
struct TemplateCond<false, If, Else>{
	typedef Else Result;
};

#endif
