/********************************************************************************
 * 
 * SDD Matrix Solver
 * 
 * Solving a SDD matrix can be reduced to solving a Laplacian matrix of size n*2.
 * This file uses the one-step reduction given in [KOSZ13], instead of the two-step Gremban reduction.
 * 
 * Public API:
 * 	SDDSolver(const Mat &A, Mat &B):
 * 		Inputs: 
 * 			Mat A: the SDD Matrix
 * 		
 * 		Outputs:
 * 			Mat B: a Laplacian matrix of size n*2
 * 			Returns a function f, which takes an AbstractSolver as input and returns an AbstractSolver
 * 			When f is fed with an AbstractSolver for B, it returns an AbstractSolver for A.
 * 
 * 	SDDSolver(const Mat &A, const function<AbstractSolver(const Mat&)> &f):
 * 		Inputs: 
 * 			Mat A: the SDD Matrix:
 * 			Function f: a function which, on input Laplacian matrix B, returns an AbstractSolver for B
 *		
 *		Outputs:
 *			An AbstractSolver for A.
 *
 *		Example:
 *			AbstractSolver S=SDDSolver(A,[](const Mat &B) {
 *				GraphSP g=IO::convertLaplacianMatrixToGraphSP(B);	//construct graph corresponding to B
 *				AbstractSolver ret=RecursiveCG(g);				//construct the Recursive CG solver
 *				B.freeMemory();							//B itself is no longer used, free memory
 *				return ret;
 *			}); 
 *
 ********************************************************************************/

#ifndef __SDD_H__
#define __SDD_H__

#include "common.h"
#include "matrix.h"
#include "graph.h"
#include "treefinder.h"
#include "abstractsolver.h"

function<AbstractSolver(const AbstractSolver&)> SDDSolver(const Mat &A, Mat &B)
{
	static FLOAT128 *d1=new FLOAT128[maxn], *d2=new FLOAT128[maxn];
	int n=A.n; 
	rep(i,0,n-1) d1[i]=0; 
	rep(i,0,n-1) d2[i]=0;
	B=Mat(n*2,n*2);
	rept(it,A.values)
	{
		int x=it->x, y=it->y; FLOAT z=it->z;
		if (x==y)
			d2[x]+=z;
		else
		{
			if (z<0)
			{
				B.entryAddValue(x,y,z);
				B.entryAddValue(x+n,y+n,z);
			}
			else if (z>0)
			{
				B.entryAddValue(x+n,y,-z);
				B.entryAddValue(x,y+n,-z);
			}
			d1[x]+=myfabs(z);
		}
	}
	rep(i,0,n-1)
	{
		FLOAT excess=d2[i]-d1[i];
		if (excess<-1e-6) 
			printf("Warning: Given matrix is not SDD, its diagonal is not dominating!\n"); 
		else if (excess<0) 
			printf("Warning: Given matrix is not SDD, its diagonal is not dominating! (floating error?)\n"); 

		B.entryAddValue(i,i,d1[i]+excess/2);
		B.entryAddValue(i+n,i+n,d1[i]+excess/2);
		B.entryAddValue(i+n,i,-excess/2);
		B.entryAddValue(i,i+n,-excess/2);
	}
	B.sortup();
	
	return [n](const AbstractSolver &S) {
		return AbstractSolver([n,S](const Vec &b, FLOAT tol, int maxit, const Vec &x0) {
			assert(n==b.n);
			Vec b1(n*2); Vec x1(n*2);
			rep(i,0,n-1) { b1[i]=b[i]; b1[i+n]=-b[i]; }
			rep(i,0,n-1) { x1[i]=x0[i]; x1[i+n]=-x0[i]; }
			SolverReturnValue ret=S.solve(b1,tol,maxit,x1);
			Vec x(n);
			rep(i,0,n-1) x[i]=(ret.x[i]-ret.x[i+n])/2;
			ret.x=x;
			return ret;
		});
	};
}

AbstractSolver SDDSolver(const Mat &A, const function<AbstractSolver(const Mat&)> &Sv) 
{
	Mat B; auto fn=SDDSolver(A,B);
	return fn(Sv(B));
}

#endif
