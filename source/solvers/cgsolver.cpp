#include "../util/io.h"
#include "../util/graph.h"
#include "cg.h"
#include "jacobiprecon.h"

int main(int argc, char *argv[])
{
	if (argc!=3) { printf("Usage: ./cgsolver graphFile rhsFile\n"); return 0; }
	
	string graphFile=argv[1];
	string rhsFile=argv[2];
	
	GraphSP g=IO::readGraphSP(graphFile);
	Mat A=IO::constructMatrixFromGraph(g);
	Vec b=IO::readMMVec(rhsFile);
	clock_t t_start = clock();
	//Vec x=CG::pcg(A,b,1e-5,-1,JacobiPreconditioner(A));
        SolverReturnValue ret=CG::pcg(A,b,1e-5,-1,JacobiPreconditioner(A));
	Vec x = ret.x;
        clock_t t_end = clock();
	double tcost=double(t_end - t_start) / double(CLOCKS_PER_SEC);
	double relres=(A*x-b).norm()/b.norm();
	//assert(relres<(1e-5)+(1e-9));
	printf("pcg %0.6e %d %.3lf\n",relres,ret.iter,tcost);
	return 0;
}
