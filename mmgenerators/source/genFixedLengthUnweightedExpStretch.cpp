#include <iostream>
#include <iomanip>
#include <math.h>

int main(int argc, char *argv[]) {

  if(argc < 4) {
    std::cerr << "Please specify n, hop, and seed." << std::endl;
    return -1;
  }
  
  size_t n = atol(argv[1]);
  size_t hop = atol(argv[2]);
 
  if(hop > n) {
    std::cerr << "hop length of paths need to be < n" << std::endl;
    return -1;
  }

  size_t precdigits=6;

  double *diag = new double[n];
  for(int i=0; i < n; ++i) {
    diag[i]=0;
  }

  double *rS = new double[n];
  
  int idx=0;

  size_t seed = atol(argv[3]);

  srand(seed);
  int m=n-hop;

  double stretch=0;

  double *printVal = new double[m];

  rS[0]=0;
  for(int i=0; i < n-1; ++i) {
    diag[i]+=1;
    diag[i+1]+=1;
    rS[i+1]=rS[i]+1;
  }

  for(int i=0; i < n-hop; ++i) {
    int u=i;
    int v=i+hop;
    double r=rS[v]-rS[u];
    while(double(rand())/double(RAND_MAX) < 0.5) {
      r/=2;
    }
    stretch+=(rS[v]-rS[u])/r;
    double roundr=round(pow(10.,precdigits)/r)/pow(10.,precdigits);
    printVal[idx]=-roundr;
    diag[u]+=roundr;
    diag[v]+=roundr;
    idx++;
  }


  printf("%%%%MatrixMarket matrix coordinate real symmetric\n%%\n%%Total Stretch %f\n", stretch);
  std::cout << n << ' ' << n  << ' ' << m+n-1+n << std::endl;


  for(int i=0; i < n-1; ++i) {
    std::cout << i+1 << ' ' << i+2 << ' ' << -1 << std::endl; 
  }

  idx=0;
  for(int i=0; i < n-hop; ++i) {
    int u=i;
    int v=i+hop;
    std::cout << u+1 << ' ' << v+1 << ' ' << printVal[idx] << std::endl;
    idx++;
  }


  for(int i=0; i < n; ++i) {
    std::cout << i+1 << ' ' << i+1 << ' ' << std::setprecision(precdigits+1) << diag[i] << std::endl;
  }

  delete diag;
  delete rS;
  delete printVal;

  return 0;
}

