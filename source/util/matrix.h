/********************************************************************************
 * 
 * Basic Definitions for Matrices and Vectors
 * 
 * Public API for Vec:
 * 	Vec(n): 			length-n vector, initialized to 0
 * 	.n:				length of vector
 *	[index]:			access the index-th coodinate (0-indexed)
 * 	.norm():			norm of vector
 *	Vec*FLOAT:			scale multiplication
 * 	Vec+Vec (Vec-Vec):	vector add/subtract
 * 	Vec*Vec:			dot multiplication
 * 
 * Public API for Mat:
 * 	Mat(n,m):			n x m sparse matrix initialized to 0
 * 	.n (.m):			size of matrix
 *	.entryAddValue(x,y,z):	Add z to entry (x,y)
 * 	.sortup():			sort the entries for better cache efficiency
 * 	.freeMemory():		destroys matrix and frees memory
 * 	Mat*Vec:			Matrix multiplication
 * 	Vec*Mat:			Matrix multiplication, Vec is treated as row vector
 * 	
 * NOTE:
 * 	Mat behaves like objects in python or javascript. That is, if you do 
 * 		Mat A; Mat B=A; 
 * 	Then A and B will be actually pointing to the same object, 
 * 	i.e. modifying B will result in A being modified as well!
 * 
 * 	However, C++ does not have an automatic garbage collection system, 
 * 	so you have to run .freeMemory() to free memory of matrices that are no longer needed.
 * 
 ********************************************************************************/

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "common.h"

struct Vec
{
	int n; FLOAT *a;
	
	Vec() { n=0; a=NULL; }
	Vec(int _n) { 
		n=_n; a=new FLOAT[n]; 
#ifdef USE_MPFR
		rep(i,0,n-1) a[i]=0.0;
#else
		memset(a,0,sizeof(FLOAT)*n); 
#endif
	}
	Vec(const Vec &b) { 
		n=b.n; a=new FLOAT[n]; 
#ifdef USE_MPFR
		rep(i,0,n-1) a[i]=b.a[i];
#else
		memcpy(a,b.a,n*sizeof(FLOAT)); 
#endif
	}
	~Vec() { n=0; delete[] a; a=NULL; }
	
	FLOAT &operator [](int i) const {
#ifndef NO_RANGE_CHECK
		if (i<0 || i>=n) assert(0);
#endif
		return a[i];
	}
	
	Vec &operator =(const Vec &b) {
		if (!a || b.n!=n)
		{
			delete[] a; n=b.n; a=new FLOAT[n];
		}
#ifdef USE_MPFR
		rep(i,0,n-1) a[i]=b.a[i];
#else
		memcpy(a,b.a,n*sizeof(FLOAT));
#endif
		return (*this);
	}
	
	FLOAT norm() const;
};

struct Mat
{
	struct MatElementType
	{
		int x,y; FLOAT z;
		bool operator <(const MatElementType &b) const
		{
			if (this->x!=b.x) return this->x<b.x;
			if (this->y!=b.y) return this->y<b.y;
			return this->z<b.z;
		}
	};
	
	int n,m;
	vector<MatElementType> &values;
	
	Mat(): n(0), m(0), values(*new vector<MatElementType>) {}
	Mat(int n, int m) : n(n), m(m), values(*new vector<MatElementType>) {}
	Mat(const Mat &b) : n(b.n), m(b.m), values(b.values) {}
	Mat& operator =(const Mat&& b) {
		n=b.n; m=b.m; values.swap(b.values);
		return (*this);
	}
	
	void entryAddValue(int x, int y, FLOAT z) {
#ifndef NO_RANGE_CHECK
		if (x<0 || x>=n || y<0 || y>=m) assert(0);
#endif
		MatElementType t; t.x=x; t.y=y; t.z=z;
		values.push_back(t);
	}
	
	void sortup() const
	{
		int all=0;
		sort(values.begin(), values.end());
		MatElementType last; last.x=-1;
		rept(it,values)
			if (it->x!=last.x || it->y!=last.y)
			{
				if (last.x>=0) { values[all]=last; all++; }
				last=(*it);
			}
			else	last.z+=it->z;
		if (last.x>=0) { values[all]=last; all++; }
		values.resize(all);
	}
	
	Mat transpose() const
	{
		Mat B(n,m);
		rept(it,values) 
		{
			MatElementType t; t.x=it->y; t.y=it->x; t.z=it->z;
			B.values.push_back(t);
		}
		B.sortup();
		return B;
	}
	
	void freeMemory() const
	{
		delete &values;
	}
};

FLOAT operator *(const Vec &a, const Vec &b)
{
	assert(a.n==b.n);
	FLOAT ret=0;
	rep(i,0,a.n-1) ret+=a.a[i]*b.a[i];
	return ret;
}

Vec operator *(const Vec &a, FLOAT b)
{
	Vec ret(a.n);
	rep(i,0,a.n-1) ret[i]=a.a[i]*b;
	return ret;
}

Vec operator +(const Vec &a, const Vec &b)
{
	assert(a.n==b.n);
	Vec ret(a.n);
	rep(i,0,a.n-1) ret[i]=a.a[i]+b.a[i];
	return ret;
}

Vec operator -(const Vec &a, const Vec &b)
{
	assert(a.n==b.n);
	Vec ret(a.n);
	rep(i,0,a.n-1) ret[i]=a.a[i]-b.a[i];
	return ret;
}

Vec operator *(const Mat &a, const Vec &b)
{
	assert(a.m==b.n);
	Vec ret(a.n);
	rept(it,a.values) ret.a[it->x]+=it->z*b.a[it->y];
	return ret;
}

Vec operator *(const Vec &a, const Mat &b)
{
	assert(a.n==b.n);
	Vec ret(b.m);
	rept(it,b.values) ret.a[it->y]+=it->z*a.a[it->x];
	return ret;
}

FLOAT Vec::norm() const
{
	return mysqrt((*this)*(*this));
} 

#endif