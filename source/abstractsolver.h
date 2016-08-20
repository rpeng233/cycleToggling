/********************************************************************************
 * 
 * Wrapper Class for a Solver
 * 
 * An "AbstractSolver" is a wrapper which hardcodes a certain matrix A, 
 * and provides the following API interface:
 * 	x = solve(b,tol,maxit,x0)
 * 	tie(x,flag) = solve(b,tol,maxit,x0)
 * 	tie(x,flag,relres) = solve(b,tol,maxit,x0)
 *	tie(x,flag,relres,iter) = solve(b,tol,maxit,x0)
 * 	tie(x,flag,relres,iter,resvec) = solve(b,tol,maxit,x0)
 * 
 * Inputs:
 * 	Vec b: the RHS of Ax=b
 * 	FLOAT tol: optional, desired tolerance, default 1e-6
 *	int maxit: optional, maximal iteration #, default -1 (unlimited)
 * 	Vec x0: optional, initial guess
 * 
 * Outputs:
 * 	Vec x: the solution vector giving smallest relative error achieved in all iterations
 * 	int flag: 
 * 		0=convergence 
 * 		1=fail to convergent to desired tolerance in maximum iteration specified
 * 		2=stagnated
 * 	FLOAT relres: the smallest relative error achieved in all iterations
 * 	int iter: actual # of iterations performed. If flag=1, this is the # of iteration in which relres is achieved
 * 	vector<FLOAT> resvec: list of relative error in every iteration
 * 
 ********************************************************************************/

#ifndef __ABSTRACTSOLVER_H__
#define __ABSTRACTSOLVER_H__

#include "common.h"

struct SolverReturnValue
{
	Vec x; int flag; FLOAT relres; int iter; vector<FLOAT> resvec;
	operator tuple<Vec&,int&,FLOAT&,int&,vector<FLOAT>& >() { return tie(x,flag,relres,iter,resvec); }
	operator tuple<Vec,int,FLOAT,int,vector<FLOAT> >() { return make_tuple(x,flag,relres,iter,resvec); }
	operator tuple<Vec&,int&,FLOAT&,int&>() { return tie(x,flag,relres,iter); }
	operator tuple<Vec,int,FLOAT,int>() { return make_tuple(x,flag,relres,iter); }
	operator tuple<Vec&,int&,FLOAT&>() { return tie(x,flag,relres); }
	operator tuple<Vec,int,FLOAT>() { return make_tuple(x,flag,relres); }
	operator tuple<Vec&,int&>() { return tie(x,flag); }
	operator tuple<Vec,int>() { return make_tuple(x,flag); }
	operator Vec() { 
		if (flag==1) fprintf(stderr,"Warning: AbstractSolver failed to converge to desired tolerance.\n"); 
		return x; 
	}
};
	
struct AbstractSolver
{
	int ty;
	function<SolverReturnValue(const Vec&,FLOAT,int,const Vec&)> fn0;
	function<SolverReturnValue(const Vec&,FLOAT,int)> fn1;
	AbstractSolver(const function<SolverReturnValue(const Vec&,FLOAT,int,const Vec&)> &fn): ty(0), fn0(fn) {};
	AbstractSolver(const function<SolverReturnValue(const Vec&,FLOAT,int)> &fn): ty(1), fn1(fn) {};
	SolverReturnValue solve(const Vec &b, FLOAT tol=1e-6, int maxit=-1, const Vec &x0=Vec(0)) const 
	{ 
		if (ty)
			return fn1(b,tol,maxit);
		else  if (x0.n==0) 
				return fn0(b,tol,maxit,Vec(b.n)); 
			else 
			{ 
				assert(x0.n==b.n); 
				return fn0(b,tol,maxit,x0); 
			}
	}
};

#endif

