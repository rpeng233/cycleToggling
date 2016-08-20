/********************************************************************************
 * 
 * Jacobi Preconditioner
 * 
 * Public API:
 * 	JacobiPreconditioner(A): Constructor.
 * 		A: the matrix
 * 
 * Example:
 * 	x = pcg(A, b, 1e-6, -1, JacobiPreconditoner(A))
 * 
 ********************************************************************************/

#ifndef __JACOBIPRECON_H__
#define __JACOBIPRECON_H__

#include "common.h"
#include "matrix.h"

struct JacobiPreconditioner
{
	int n;
	Vec diag;
	
	JacobiPreconditioner(const Mat &A)
	{
		n=A.n;
		Vec d(n);
		rept(it,A.values) if (it->x==it->y) d[it->x]+=it->z;
		diag=d;
	}
	
	Vec operator()(const Vec &b)
	{
		assert(b.n==n);
		Vec ret(n);
		rep(i,0,n-1) ret[i]=b[i]/diag[i];
		return ret;
	}
};

#endif
