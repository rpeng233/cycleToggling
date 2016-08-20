#ifndef __RICHARDSON_H__
#define __RICHARDSON_H__

#include "common.h"
#include "matrix.h"
#include "abstractsolver.h"


AbstractSolver Richardson(const Mat &A, const GraphSP &g)
{
	AbstractSolver S=KOSZ(UltraSparsifier(g,400));
	return AbstractSolver([A,S](const Vec &b, FLOAT tol, int maxit, const Vec &x0) {
		Vec x=x0;
		FLOAT bnorm=b.norm();
		while (1)
		{
			Vec r=b-A*x;
			FLOAT relres=r.norm()/bnorm;
			printf("%.16lf\n",relres);
			if (relres<tol) { SolverReturnValue ret; ret.x=x; ret.flag=0; return ret; }
			x=x+S.solve(r)*FLOAT(0.01);
		}
	});
}

AbstractSolver Richardson(const Mat &A, const AbstractSolver &precon)
{
	return AbstractSolver([A,precon](const Vec &b, FLOAT tol, int maxit, const Vec &x0)->SolverReturnValue {
		Vec x=x0; Vec bestx; FLOAT bestrelres; int whichit;
		Vec r0=precon.solve(b-A*x);
		int i=0; FLOAT bnorm=b.norm();
		SolverReturnValue ret;
		vector<FLOAT> &rm=ret.resvec;
		while (1)
		{
			Vec ar0=A*r0;
			Vec mar0=precon.solve(ar0);
			FLOAT nr0=r0*r0;
			FLOAT relres=mysqrt(nr0)/bnorm;
			//rm.push_back(relres);
			//printf("Richardson Iteration #%d, relres %.16lf\n",i,relres);
			if (relres<tol) { ret.x=x; ret.flag=0; ret.iter=i; ret.relres=relres; return ret; }
			if (maxit!=-1 && (i==0 || relres<bestrelres)) { bestrelres=relres; bestx=x; whichit=i; }
			i++;
			if (maxit!=-1 && i>maxit) { ret.x=bestx; ret.flag=1; ret.iter=whichit; ret.relres=bestrelres; return ret; }
			FLOAT beta=r0*mar0;
			FLOAT alpha=nr0/beta;
			x=x+(r0*alpha);
			//r0=r0-(mar0*alpha);
			r0=precon.solve(b-A*x);
		}
	});
}

#endif