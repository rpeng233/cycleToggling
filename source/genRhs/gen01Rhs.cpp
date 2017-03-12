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
	Vec b(n);
        b[0]=-1;
        b[n-1]=1;
	rep(i,1,n-2) b[i]=0;
	
	//string gradfname=fname.substr(0,fname.length()-4)+"_01_rhs.mtx";
	IO::saveMMVec(b,outname);
}

int main(int argc, char *argv[])
{
	if (argc<3) { printf("needs filename\n"); return 0; }
	lemon(argv[1],argv[2]);
	return 0;
}

