#include <cstdlib>
#include <cstdio> 
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;

double diag[1000000];

int main(int argc, char *argv[])
{
	if (argc < 2) {
		cerr << "Please specify n." << endl;
		return -1;
	}

	size_t n = atol(argv[1]);
	
	int m = n-2;

	printf("%%%%MatrixMarket matrix coordinate real symmetric\n%%\n");
	cout << n << ' ' << n  << ' ' << m+n+n << endl;


	for (size_t i = 0; i < n - 1; i++) {
		cout << i+1 << ' ' << i + 2 << ' ' << -1 << endl;
		diag[i]+=1;
		diag[i+1]+=1;
	}

	size_t s = 0;
	for (size_t i = 0; i < m; i++) {
		double r=2;
		while (double(rand())/double(RAND_MAX)<0.5) r/=2;
		cout << i+1 << ' ' << i+3 << ' ' << -1./r << endl;
		diag[i]+=1./r;
		diag[i+2]+=1./r;
	}

	for(size_t i=0; i < n; ++i) {
	  cout << i+1 << ' ' << i+1 << ' ' << diag[i] << std::endl;
	}

	return 0;
}
