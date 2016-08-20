//#define BENCHMARK_MODE
#include "common.h"
#include "matrix.h"
#include "recursivecg.h"
#include "io.h"
#include "sddmsolver.h"
#include "jacobiprecon.h"
#include "treeprecon.h"


void lemon(string fname)
{
	GraphSP g=IO::readGraphSP(fname);
	Mat A=IO::constructMatrixFromGraph(g);
	int n=A.n;
	Vec b(n);
        b[0]=-1;
        b[n-1]=1;
	rep(i,1,n-2) b[i]=0;
	
	string gradfname=fname.substr(0,fname.length()-4)+"_01_rhs.mtx";
	IO::saveMMVec(b,gradfname);
}

int main(int argc, char *argv[])
{
	if (argc<2) { printf("needs filename\n"); return 0; }
	lemon(argv[1]);
	return 0;
}

