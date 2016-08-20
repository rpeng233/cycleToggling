#ifndef __CHOL_H__
#define __CHOL_H__

#include "common.h"
#include "matlab.h"
#include "matrix.h"
#include "io.h"
#include "abstractsolver.h"

namespace CholHelper
{
	//solve R'*x=b
	//R' is upper triangular (so R itself should be lower triangular!)
	//This assumes R is sorted!
	Vec upper_back_substitution(const Mat &R, const Vec &_b)
	{
		assert(R.n==_b.n);
		Vec ret(R.n), b=_b;
		for (auto it=R.values.rbegin(); it!=R.values.rend(); it++)
		{
			int x=it->x, y=it->y; FLOAT z=it->z;
			if (x==y) ret.a[x]=b.a[x]/z; else b.a[y]-=ret.a[x]*z;
		}
		return ret;
	}
	
	//solve R'*x=b
	//R' is lower triangular (so R itself should be upper triangular!)
	//This assumes R is sorted!
	Vec lower_back_substitution(const Mat &R, const Vec &_b)
	{
		assert(R.n==_b.n);
		Vec ret(R.n), b=_b;
		rept(it,R.values)
		{
			int x=it->x, y=it->y; FLOAT z=it->z;
			if (x==y) ret.a[x]=b.a[x]/z; else b.a[y]-=ret.a[x]*z;
		}
		return ret;
	}
}

AbstractSolver CholLaplacianSolver(const Mat &A, const function<int(int)> decision = NULL)
{
	string fA=Matlab::gettmpfile();
	string fR=Matlab::gettmpfile();
	string fS=Matlab::gettmpfile();
	IO::saveMM(A,fA);
	string cmd="tic;A=mmread('"+fA+"');toc; n=size(A,1); A=A(1:n-1,1:n-1); tic;[R,p,S]=chol(A);toc; tic;mmwrite('"+fR+"',R);toc; tic;mmwrite('"+fS+"',S);toc;";
	//string cmd="tic;A=mmread('"+fA+"');toc; n=size(A,1); A=A(1:n-1,1:n-1); tic;R=ichol(A, struct('type','ict','droptol',1e-2));toc; S=speye(n-1); tic;mmwrite('"+fR+"',R);toc; tic;mmwrite('"+fS+"',S);toc;";
	Matlab::execute(cmd);
	Mat R=IO::readMMnonsym(fR);
	Mat S=IO::readMMnonsym(fS);
	assert(R.m==R.n); assert(R.n==A.n-1);
	string rmcmd="rm "+fA+" "+fR+" "+fS;
	system(rmcmd.c_str());
	
	if (decision!=NULL && !(decision(R.values.size())))
	{
		R.freeMemory(); S.freeMemory();
		return AbstractSolver([](const Vec &_b, double tol, int iter)->SolverReturnValue { assert(0); });
	}
	/*
	 * R'*R=S'*A*S
	 * Ax=b x=Sy
	 * S'A(Sy)=S'b
	 * R'Ry=S'b
	 */
	Mat RT=R.transpose();
	Mat ST=S.transpose();
	return AbstractSolver([A,R,RT,S,ST](const Vec &_b, double tol, int iter) {
		assert(_b.n==R.n+1);
		Vec b(_b.n-1); rep(i,0,b.n-1) b.a[i]=_b.a[i];
		b=ST*b;
		Vec x=S*CholHelper::upper_back_substitution(RT,CholHelper::lower_back_substitution(R,b));
		SolverReturnValue ret; ret.x=Vec(b.n+1); 
		rep(i,0,b.n-1) ret.x.a[i]=x.a[i];
		ret.flag=0; ret.relres=0;
		printf("Chol relres: %.16lf\n",(A*ret.x-_b).norm()/_b.norm());	
		return ret;
	});
}

#endif
