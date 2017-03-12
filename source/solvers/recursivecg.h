/********************************************************************************
 * 
 * Recursive CG Graph Laplacian Solver
 * 
 * Public API:
 * 	RecursiveCG(const GraphSP &g, int flag): 
 * 		Inputs:
 * 			GraphSP g: the graph
 * 			int flag: the RecursiveCG solver will call JacobiCG with max=flag iters, 
 * 			          to get a good initial solution (default 10)
 * 	
 * 		Returns:
 * 			an AbstractSolver for g
 * 			see abstractsolver.h for definition
 * 
 * FIXME: 
 * 	atm the shrinkStrategy is very crude, can tweak for better performance
 * 
 ********************************************************************************/

#ifndef __RECURSIVECG_H__
#define __RECURSIVECG_H__

#include "../util/common.h"
#include "../util/matrix.h"
#include "../util/graph.h"
#include "cg.h"
#include "kosz.h"
#include "../util/io.h"
#include "ultrasparsifier.h"
#include "contractsolver.h"
#include "jacobiprecon.h"
#include "treeprecon.h"
#include "../util/treefinder.h"
#include "richardson.h"
#include "chol.h"


namespace RecursiveCGHelper
{
	int shrinkStrategy(const GraphSP &g)	//returns # of OTE to be sampled, -1 for direct solver
	{
		int n=g.n, sm;
		if (n>=10000000) 
			sm=50000;
		else if (n>=100000)
			sm=5000;
		else sm=-1;
		return sm;
	}
}

AbstractSolver RecursiveCG(const GraphSP &g, int flag = 10)
{
	int sm=RecursiveCGHelper::shrinkStrategy(g);
	//if (sm==-1) return Richardson(IO::constructMatrixFromGraph(g),g);
	//if (sm==-1) { Mat mG=IO::constructMatrixFromGraph(g); return PCG(mG,TreePreconditioner(g)); }
	if (sm==-1) { 
		int flag = 1;
		return CholLaplacianSolver(IO::constructMatrixFromGraph(g), [&flag](int sz) mutable -> int {
#ifdef BENCHMARK_MODE
			sprintf(mysql_buffer,"cholnnz = '%d'",sz);
			Mysql::updaterow();
#endif		
			return 1;
		}); 
	}
	//GraphSP g2=UltraSparsifier(g,sm); 
	//GraphSP g2=UltraSparsifier(g,1,(g.o.size()+g.n)/20);
	GraphSP g2=UltraSparsifier(g,0.01,(g.n==1000000?g.n/20:g.n/40));
	AbstractSolver preconS=ContractSolver(g2,[flag](const GraphSP &h) { 
#ifdef BENCHMARK_MODE
		sprintf(mysql_buffer,"chainsize = CONCAT(chainsize,'%d,'), chainote = CONCAT(chainote,'%d,')",h.n,h.o.size());
		Mysql::updaterow();
#endif
		return RecursiveCG(h,0); 
	});
	Mat mG=IO::constructMatrixFromGraph(g);
	//Mat mG2=IO::constructMatrixFromGraph(g2);
	g2.freeMemory();
	//AbstractSolver S1=PCG(mG,preconS);
	
	AbstractSolver S1=PCG(mG,[preconS/*,mG2*/](const Vec &b, const Vec &x, FLOAT err, FLOAT tol) mutable {
		Vec xs; int flag; 
		tie(xs,flag)=preconS.solve(b,0.5,50);	//solve to 0.1 tolerance only
		//printf("Actual relres: %.16lf\n",(mG2*xs-b).norm()/b.norm());
		return xs;
	});
	
	if (flag>0)
	{
		int ff=flag;
		AbstractSolver S=PCG(mG,JacobiPreconditioner(mG));
		//AbstractSolver S=PCG(mG,TreePreconditioner(g));
		return AbstractSolver([S,S1,ff](const Vec &b, FLOAT tol, int maxit, const Vec &x0) {
			//iterate ff times for a good starting solution
			Vec x1; int flag; tie(x1,flag)=S.solve(b,tol,ff,x0);
			return S1.solve(b,tol,maxit,x1);
		});
	}
	else	return S1;
}

#endif
