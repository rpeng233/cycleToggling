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
	Vec x(n);
	rep(i,0,n-1) x[i]=(double(rand())/double(RAND_MAX)-0.5)*10;
	Vec b=A*x;
	string gradfname=fname.substr(0,fname.length()-4)+"_Rand_rhs.mtx";
	IO::saveMMVec(b,gradfname);
}

int main(int argc, char *argv[])
{
	if (argc<2) { printf("needs filename\n"); return 0; }
	lemon(argv[1]);
	return 0;
}

