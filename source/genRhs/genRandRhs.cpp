//#define BENCHMARK_MODE
#include "../util/common.h"
#include "../util/matrix.h"
#include "../solvers/recursivecg.h"
#include "../util/io.h"
#include "../solvers/sddmsolver.h"
#include "../solvers/jacobiprecon.h"
#include "../solvers/treeprecon.h"


void lemon(string fname, string outname)
{
	GraphSP g=IO::readGraphSP(fname);
	Mat A=IO::constructMatrixFromGraph(g);
	int n=A.n;
	Vec x(n);
	rep(i,0,n-1) x[i]=(double(rand())/double(RAND_MAX)-0.5)*10;
	Vec b=A*x;
	//string gradfname=fname.substr(0,fname.length()-4)+"_Rand_rhs.mtx";
	IO::saveMMVec(b,outname);
}

int main(int argc, char *argv[])
{
	if (argc<3) { printf("needs filename\n"); return 0; }
	lemon(argv[1],argv[2]);
	return 0;
}

