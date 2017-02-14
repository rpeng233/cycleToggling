#include <iostream>
#include <iomanip>
#include <math.h>

int main(int argc, char *argv[])
{
  if (argc < 2) {
    std::cerr << "Please specify n." << std::endl;
    return -1;
  }

  size_t n = atol(argv[1]);
  
  size_t precdigits=6;

  double *diag = new double[n];
  for(int i=0; i < n; ++i) {
    diag[i]=0;
  }
	
  int m=n-2;

  printf("%%%%MatrixMarket matrix coordinate real symmetric\n%%\n");
  std::cout << n << ' ' << n  << ' ' << m+n-1+n << std::endl;


  for(int i=0; i < n-1; ++i) {
    std::cout << i+1 << ' ' << i + 2 << ' ' << -1 << std::endl;
    diag[i]+=1;
    diag[i+1]+=1;
  }

  for(int i=0; i < m; ++i) {
    double r=2;
    while(double(rand())/double(RAND_MAX) < 0.5) {
      r/=2;
    }
    double roundr=round(pow(10.,precdigits)/r)/pow(10.,precdigits);
    std::cout << i+1 << ' ' << i+3 << ' ' << -roundr << std::endl;
    diag[i]+=roundr;
    diag[i+2]+=roundr;
  }
  
  for(int i=0; i < n; ++i) {
    std::cout << i+1 << ' ' << i+1 << ' ' << std::setprecision(precdigits+1) << diag[i] << std::endl;
  }

  delete diag;
  return 0;
}
