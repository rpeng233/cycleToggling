#ifndef __COMMON_H__
#define __COMMON_H__

//uncomment the following macro to use __float128 instead of double
//remember to add -lquadmath to compile options
//#define USE_FLOAT128	

//uncomment the following macro to use MPFR instead of double
//remember to add -Wall -ansi -pedantic -lmpfr to compile options
//overrides USE_FLOAT128
//#define USE_MPFR
#define MPFR_PRECISION 1024	//bit precision of MPFR

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <cstring>
#include <tuple>
#include <cassert>
#include <functional>
#include <iomanip> 
#include <fstream>
#include <sstream>
#ifdef USE_MPFR
#include "mpfrreal.hpp"
#else
#ifdef USE_FLOAT128
extern "C" {
	#include <quadmath.h>
}
#endif
#endif

using namespace std;

typedef long long LL;
typedef unsigned long long ULL;

#define SIZE(x) (int((x).size()))
#define rep(i,l,r) for (int i=(l); i<=(r); i++)
#define repd(i,r,l) for (int i=(r); i>=(l); i--)
#define rept(i,c) for (__typeof((c).begin()) i=(c).begin(); i!=(c).end(); i++)

#ifndef ONLINE_JUDGE
#define debug(x) { cerr<<#x<<" = "<<(x)<<endl; }
#else
#define debug(x) {}
#endif

#define maxn 1050001

#ifdef USE_MPFR
typedef mpfr::real<MPFR_PRECISION> FLOAT;
typedef mpfr::real<MPFR_PRECISION> FLOAT128;
#define mysqrt sqrt
#define myfabs fabs
double printFloat(FLOAT x) { 
	stringstream s; s<<setprecision(20)<<scientific<<x;
	double ret; sscanf(s.str().c_str(),"%lf",&ret);
	return ret; 
}
#else
#ifdef USE_FLOAT128
typedef __float128 FLOAT;
typedef __float128 FLOAT128;
#define mysqrt sqrtq
#define myfabs fabsq
double printFloat(FLOAT x) { return (double)x; }
#else
typedef double FLOAT;
typedef __float128 FLOAT128;
#define mysqrt sqrt
#define myfabs fabs
double printFloat(FLOAT x) { return x; }
#endif
#endif

#define formatf(len,prec,value) setw(len)<<setprecision(prec)<<fixed<<value
#define formate(len,prec,value) setw(len)<<setprecision(prec)<<scientific<<value
#define formatd(len,value) setw(len)<<fixed<<value
#define formats(len,value) setw(len)<<value

#endif
