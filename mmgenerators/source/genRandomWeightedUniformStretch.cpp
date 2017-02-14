#include <iostream>
#include <iomanip>
#include <math.h>
#include <set>

std::set<std::pair<int,int> > has;

int main(int argc, char *argv[]) {

  if(argc < 3) {
    std::cerr << "Please specify n and m" << std::endl;
  }
  size_t n=atol(argv[1]);
  size_t m=atol(argv[2]);
  size_t precdigits=6;
  if(argc == 4) {
    precdigits=atol(argv[3]);
  }

  double *rS = new double[n];

  double *diag = new double[n];
  for(int i=0; i < n; ++i) {
    diag[i]=0;
  }

  printf("%%%%MatrixMarket matrix coordinate real symmetric\n%%\n%%Total Stretch %d\n", (int)m);
  std::cout << n << ' ' << n  << ' ' << m+n-1+n << std::endl;

  rS[0]=0;
  for(int i=0; i < n-1; ++i) {
    double r=rand()%1000+1;
    double roundr=round(pow(10.,precdigits)/r)/pow(10.,precdigits);
    std::cout << i+1 << ' ' << i+2 << ' ' << -roundr << std::endl;
    diag[i]+=roundr;
    diag[i+1]+=roundr;
    rS[i+1]=rS[i]+r;
  }

  has.clear();

  for(int i=0; i < m; ++i) {
    int u=rand()%n;
    int v=rand()%n;

    if(u > v) {
      std::swap(u, v);
    }

    while(u == v || v == u+1 || has.find(std::make_pair(u, v)) != has.end()) {
      u=rand()%n;
      v=rand()%n;
      if(u > v) {
        std::swap(u, v);
      }
    }

    has.insert(std::make_pair(u, v));

   
    double roundr=round(pow(10.,precdigits)/(rS[v]-rS[u]))/pow(10.,precdigits);
    if(roundr == 0) {
      std::cerr << "increase precision because edge weights too small" << std::endl;
      return -1;
    }
    std::cout << u+1 << ' ' << v+1 << ' ' << -roundr << std::endl;

    diag[u]+=roundr;
    diag[v]+=roundr;
  }


  for(int i=0; i < n; ++i) {
    std::cout << i+1 << ' ' << i+1 << ' ' << std::setprecision(precdigits+1) << diag[i] << std::endl;
  }

  delete diag;
  delete rS;

  return 0;
}

