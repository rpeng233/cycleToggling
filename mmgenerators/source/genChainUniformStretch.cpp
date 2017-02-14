#include <iostream>


int main(int argc, char *argv[])
{
  if (argc < 2) {
    std::cerr << "Please specify n." << std::endl;
    return -1;
  }
  
  size_t n = atol(argv[1]);

  double *diag = new double[n];
  for(int i=0; i < n; ++i) {
    diag[i]=0;
  }
        
  int m=n-2;

  printf("%%%%MatrixMarket matrix coordinate real symmetric\n%%\n%%Total Stretch %d\n", m);
  std::cout << n << ' ' << n  << ' ' << m+n-1+n << std::endl;
  
  
  for(int i=0; i < n-1; ++i) {
    std::cout << i+1 << ' ' << i+2 << ' ' << -1 << std::endl;
    diag[i]+=1;
    diag[i+1]+=1;
  }
  
  for(int i=0; i < m; ++i) {
    std::cout << i+1 << ' ' << i+3 << ' ' << -.5 << std::endl;
    diag[i]+=.5;
    diag[i+2]+=.5;
  }
  
  for(int i=0; i < n; ++i) {
    std::cout << i+1 << ' ' << i+1 << ' ' << diag[i] << std::endl;
  }
  
  delete diag;

  return 0;
}
