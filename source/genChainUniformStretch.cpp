#include <cstdlib>
#include <cstdio> 
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char *argv[])
{
	if (argc < 3) {
		cerr << "Please specify n and k." << endl;
		return -1;
	}

	size_t n = atol(argv[1]);
	size_t k = atol(argv[2]);
	
	int m = n-2;

	cout << n << ' ' << m + n - 1 << ' ' << k << endl;

	for (size_t i = 0; i < n - 1; i++) {
		cout << i << ' ' << i + 1 << ' ' << 1 << endl;
	}

	size_t s = 0;
	for (size_t i = 0; i < m; i++) {
		cout << i << ' ' << i+2 << ' ' << 2 << endl;
	}

	for (size_t i = 0; i < n - 1 + m; i++) {
		if (i < n - 1) {
			cout << "1.0" << endl;
		} else {
			cout << "0.0" << endl;
		}
	}
	/*
	for (size_t i = 0; i < k; ++i) {
		printf("%d\n", n + rand() % m);
		}*/

	return 0;
}
