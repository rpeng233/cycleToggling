#include <cstdlib>
#include <cstdio> 
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include "io.h"
#include "graph.h"
#include "matrix.h"

using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char *argv[])
{
	if (argc < 3) {
		cerr << "Please specify graph file and output file" << endl;
		return -1;
	}

	string graphFile=argv[1];
	string outFile=argv[2];
	GraphSP g=IO::readGraphSP(graphFile);
	Mat A=IO::constructMatrixFromGraph(g);
	IO::saveMM(A, outFile);
	return 0;
}
