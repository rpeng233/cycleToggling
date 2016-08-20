/********************************************************************************
 * 
 * Preconditioned Conjugacy Gradiant Method
 * 
 * Public API (almost identical as in MATLAB):
 * 	CG::pcg:
 * 		... = CG::pcg(A,b)
 * 		... = CG::pcg(A,b,tol)
 * 		... = CG::pcg(A,b,tol,maxit)
 * 		... = CG::pcg(A,b,tol,maxit,Preconditioner)
 * 		... = CG::pcg(A,b,tol,maxit,Preconditioner,x0)
 * 		LHS accepts any formats defined in abstractsolver.h
 * 
 * 	PCG(A,precon):
 * 		construct AbstractSolver
 * 		precon is optional (see below)
 * 
 * Inputs:
 * 	Mat A, Vec b: Ax=b, A square real symmetric positive definite
 * 	FLOAT tol: desired tolerance (default 1e-6)
 * 	int maxit: max # of iterations (default -1, unlimited)
 * 	Preconditioner (default none):
 * 		Either
 * 			an "AbstractSolver" type variable, 
 * 		or an function of any of the 4 prototypes:
 * 			Vec Preconditioner(const Vec &b);
 * 			Vec Preconditioner(const Vec &b, FLOAT err);
 * 			Vec Preconditioner(const Vec &b, const Vec &x, FLOAT err);
 *			Vec Preconditioner(const Vec &b, const Vec &x, FLOAT err, FLOAT tol);
 * 		Vec x is the current solution. FLOAT err is the current relative error. 
 * 		It should return M^{-1}b. 
 * 	Vec x0: initial guess (default zero vector)
 *	 
 ********************************************************************************/

#ifndef __CG_H__
#define __CG_H__

#include "common.h"
#include "matrix.h"
#include "abstractsolver.h"

namespace CG
{
	void pcgsolve(const Mat &A, const Vec &b, FLOAT tol, int maxit, const function<void(const Vec&, const Vec&, FLOAT, FLOAT, Vec&)> &precon, const Vec &_x0, SolverReturnValue &ret)
	{
		assert(A.n==A.m); assert(A.n==b.n);
		Vec &x0=ret.x; x0=_x0; Vec r0=b-A*x0;
		FLOAT bnorm=b.norm(), err=r0.norm()/bnorm; 
		vector<FLOAT> &rm=ret.resvec; rm.push_back(err);
		if (err<tol) { ret.flag=0; ret.relres=err; ret.iter=0; return; }
		Vec mr; precon(r0,x0,err,tol,mr);
		Vec d0=mr;
		FLOAT minerr=err; Vec bestx=x0; int whichit=0;
		int i=0; Vec tmp; FLOAT t=(r0*mr); FLOAT lasterr=err; int streak=0;
		while (1)
		{
			i++; if (maxit!=-1 && i>maxit) break;
			tmp=A*d0;
			FLOAT alpha=t/(d0*tmp);
			x0=x0+d0*alpha;
			r0=r0-tmp*alpha;
			FLOAT err=r0.norm()/bnorm; rm.push_back(err);
			precon(r0,x0,err,tol,mr);
			FLOAT nt=r0*mr;
			FLOAT beta=nt/t; t=nt;
			d0=mr+d0*beta;
			//printf("CG Iteration #%d, relres %.16lf\n",i,printFloat(err));
			if (maxit!=-1 && err<minerr) { minerr=err; bestx=x0; whichit=i; }
			if (err<tol) { ret.flag=0; ret.relres=err; ret.iter=i; return; }
			if (myfabs((err-lasterr)/lasterr)<1e-6)
			{
				streak++;
				if (streak>=10) { ret.flag=2; ret.relres=err; ret.iter=i; return; }
			}
			else	streak=0;
			lasterr=err;
		}
		ret.x=bestx; ret.flag=1; ret.relres=minerr; ret.iter=whichit; return;
	}

	SolverReturnValue pcg(const Mat &A, const Vec &b, FLOAT tol, int maxit, const AbstractSolver &precon, const Vec &x0)
	{
		SolverReturnValue ret;
		pcgsolve(A,b,tol,maxit,[&precon](const Vec &b, const Vec &x, FLOAT err, FLOAT tol, Vec &ret) { ret=precon.solve(b,tol); },x0,ret);
		return ret;
	}
	
	SolverReturnValue pcg(const Mat &A, const Vec &b, FLOAT tol, int maxit, const function<Vec(const Vec&, const Vec&, FLOAT, FLOAT)> &precon, const Vec &x0)
	{
		SolverReturnValue ret;
		pcgsolve(A,b,tol,maxit,[&precon](const Vec &b, const Vec &x, FLOAT err, FLOAT tol, Vec &ret) { ret=precon(b,x,err,tol); },x0,ret);
		return ret;
	}
	
	SolverReturnValue pcg(const Mat &A, const Vec &b, FLOAT tol, int maxit, const function<Vec(const Vec&, const Vec&, FLOAT)> &precon, const Vec &x0)
	{
		SolverReturnValue ret;
		pcgsolve(A,b,tol,maxit,[&precon](const Vec &b, const Vec &x, FLOAT err, FLOAT tol, Vec &ret) { ret=precon(b,x,err); },x0,ret);
		return ret;
	}
	
	SolverReturnValue pcg(const Mat &A, const Vec &b, FLOAT tol, int maxit, const function<Vec(const Vec&, FLOAT)> &precon, const Vec &x0)
	{
		SolverReturnValue ret;
		pcgsolve(A,b,tol,maxit,[&precon](const Vec &b, const Vec &x, FLOAT err, FLOAT tol, Vec &ret) { ret=precon(b,err); },x0,ret);
		return ret;
	}
	
	SolverReturnValue pcg(const Mat &A, const Vec &b, FLOAT tol, int maxit, const function<Vec(const Vec&)> &precon, const Vec &x0)
	{
		SolverReturnValue ret;
		pcgsolve(A,b,tol,maxit,[&precon](const Vec &b, const Vec &x, FLOAT err, FLOAT tol, Vec &ret) { ret=precon(b); },x0,ret);
		return ret;
	}
	
	SolverReturnValue pcg(const Mat &A, const Vec &b, FLOAT tol, int maxit, const AbstractSolver &precon)
	{
		SolverReturnValue ret;
		pcgsolve(A,b,tol,maxit,[&precon](const Vec &b, const Vec &x, FLOAT err, FLOAT tol, Vec &ret) { ret=precon.solve(b,tol); },Vec(A.m),ret);
		return ret;
	}
	
	SolverReturnValue pcg(const Mat &A, const Vec &b, FLOAT tol, int maxit, const function<Vec(const Vec&, const Vec&, FLOAT, FLOAT)> &precon)
	{
		SolverReturnValue ret;
		pcgsolve(A,b,tol,maxit,[&precon](const Vec &b, const Vec &x, FLOAT err, FLOAT tol, Vec &ret) { ret=precon(b,x,err,tol); },Vec(A.m),ret);
		return ret;
	}
	
	SolverReturnValue pcg(const Mat &A, const Vec &b, FLOAT tol, int maxit, const function<Vec(const Vec&, const Vec&, FLOAT)> &precon)
	{
		SolverReturnValue ret;
		pcgsolve(A,b,tol,maxit,[&precon](const Vec &b, const Vec &x, FLOAT err, FLOAT tol, Vec &ret) { ret=precon(b,x,err); },Vec(A.m),ret);
		return ret;
	}
	
	SolverReturnValue pcg(const Mat &A, const Vec &b, FLOAT tol, int maxit, const function<Vec(const Vec&, FLOAT)> &precon)
	{
		SolverReturnValue ret;
		pcgsolve(A,b,tol,maxit,[&precon](const Vec &b, const Vec &x, FLOAT err, FLOAT tol, Vec &ret) { ret=precon(b,err); },Vec(A.m),ret);
		return ret;
	}
	
	SolverReturnValue pcg(const Mat &A, const Vec &b, FLOAT tol, int maxit, const function<Vec(const Vec&)> &precon)
	{
		SolverReturnValue ret;
		pcgsolve(A,b,tol,maxit,[&precon](const Vec &b, const Vec &x, FLOAT err, FLOAT tol, Vec &ret) { ret=precon(b); },Vec(A.m),ret);
		return ret;
	}

	SolverReturnValue pcg(const Mat &A, const Vec &b, FLOAT tol, int maxit)
	{
		SolverReturnValue ret;
		pcgsolve(A,b,tol,maxit,[](const Vec &b, const Vec &x, FLOAT err, FLOAT tol, Vec &ret) { ret=b; },Vec(A.m),ret);
		return ret;
	}

	SolverReturnValue pcg(const Mat &A, const Vec &b, FLOAT tol)
	{
		SolverReturnValue ret;
		pcgsolve(A,b,tol,-1,[](const Vec &b, const Vec &x, FLOAT err, FLOAT tol, Vec &ret) { ret=b; },Vec(A.m),ret);
		return ret;
	}

	SolverReturnValue pcg(const Mat &A, const Vec &b)
	{
		SolverReturnValue ret;
		pcgsolve(A,b,1e-6,-1,[](const Vec &b, const Vec &x, FLOAT err, FLOAT tol, Vec &ret) { ret=b; },Vec(A.m),ret);
		return ret;
	} 
}

AbstractSolver PCG(const Mat &A, const AbstractSolver &precon)
{
	return AbstractSolver([A,precon](const Vec &b, FLOAT tol, int maxit, const Vec &x0) {
		return CG::pcg(A,b,tol,maxit,precon,x0);
	});
}

AbstractSolver PCG(const Mat &A, const function<Vec(const Vec&, const Vec&, FLOAT, FLOAT)> &precon)
{
	return AbstractSolver([A,precon](const Vec &b, FLOAT tol, int maxit, const Vec &x0) {
		return CG::pcg(A,b,tol,maxit,precon,x0);
	});
}

AbstractSolver PCG(const Mat &A, const function<Vec(const Vec&, const Vec&, FLOAT)> &precon)
{
	return AbstractSolver([A,precon](const Vec &b, FLOAT tol, int maxit, const Vec &x0) {
		return CG::pcg(A,b,tol,maxit,precon,x0);
	});
}

AbstractSolver PCG(const Mat &A, const function<Vec(const Vec&, FLOAT)> &precon)
{
	return AbstractSolver([A,precon](const Vec &b, FLOAT tol, int maxit, const Vec &x0) {
		return CG::pcg(A,b,tol,maxit,precon,x0);
	});
}

AbstractSolver PCG(const Mat &A, const function<Vec(const Vec&)> &precon)
{
	return AbstractSolver([A,precon](const Vec &b, FLOAT tol, int maxit, const Vec &x0) {
		return CG::pcg(A,b,tol,maxit,precon,x0);
	});
}

AbstractSolver PCG(const Mat &A)
{
	return AbstractSolver([A](const Vec &b, FLOAT tol, int maxit, const Vec &x0) {
		return CG::pcg(A,b,tol,maxit,[](const Vec &z) { return z; },x0);
	});
}
		
#endif
