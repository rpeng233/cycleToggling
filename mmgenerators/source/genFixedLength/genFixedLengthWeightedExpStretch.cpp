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

  if(argc == 5) {
    precdigits = atol(argv[4]);
  }

  double *diag = new double[n];
  for(int i=0; i < n; ++i) {
    diag[i]=0;
  }

  double *rS = new double[n];
  

  size_t seed = atol(argv[3]);

  srand(seed);
  int m=n-hop;

  double *printVal = new double[m];
  int idx=0;

  double stretch=0;
  
  rS[0]=0;
  for(int i=0; i < n-1; ++i) {
    double r=rand()%1000+1;
    rS[i+1]=rS[i]+r;
  }

  for(int i=0; i < n-hop; ++i) {
    int u=i;
    int v=i+hop;
    double r=rS[v]-rS[u];
    double newr=r;
    double maxr=0;
    for(int j=0; j < hop; ++j) {
      if(rS[u+j+1]-rS[u+j] > maxr) {
        maxr=rS[u+j+1]-rS[u+j];
      }
    }
    while(double(rand())/double(RAND_MAX) < 0.5) {
      if((newr/2+r-maxr)/maxr < 2*r/newr) {
        break;
      }
      newr/=2;
    }
    stretch+=(r/newr);
    if((newr+r-maxr)/maxr < r/newr) {
      std::cerr << "low stretch tree no longer path" << std::endl;
      return -1;
    }
    double roundr=round(pow(10.,precdigits)/newr)/pow(10.,precdigits);
    if(roundr == 0) {
      std::cerr << "increase precision because edge weights too small" << std::endl;
      return -1;
    }
    printVal[idx]=-roundr;
    diag[u]+=roundr;
    diag[v]+=roundr;
    idx++;
  }


  printf("%%%%MatrixMarket matrix coordinate real symmetric\n%%\n%%Total Stretch %f\n", stretch);
  std::cout << n << ' ' << n  << ' ' << m+n-1+n << std::endl;

  for(int i=0; i < n-1; ++i) {
    double roundr=round(pow(10.,precdigits)/(rS[i+1]-rS[i]))/pow(10.,precdigits);
    if(roundr == 0) {
      std::cerr << "increase precision because edge weights too small" << std::endl;
      return -1;
    }
    std::cout << i+1 << ' ' << i+2 << ' ' << -roundr << std::endl; 
    diag[i]+=roundr;
    diag[i+1]+=roundr;
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

