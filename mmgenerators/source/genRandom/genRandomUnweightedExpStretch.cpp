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

  double *printVal = new double[m];
  int idx=0;

  
  rS[0]=0;
  for(int i=0; i < n-1; ++i) {
    rS[i+1]=rS[i]+1;
  }

  has.clear();

  double stretch=0;

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
  }

  std::set<std::pair<int,int> >::iterator iter;
  for(iter=has.begin(); iter != has.end(); ++iter) {
    int u=(*iter).first;
    int v=(*iter).second;

    double r=rS[v]-rS[u];
    double newr=r;

    while (double(rand())/double(RAND_MAX) < 0.5) {
      if(newr/2+r-1 < 2*r/newr) {
        break;
      }
      newr/=2;
    }
    stretch+=(r/newr);
    if((newr+r-1) < r/newr) {
      std::cerr << "low stretch tree no longer path" << std::endl;
      return -1;
    }
    double roundr=round(pow(10.,precdigits)/newr)/pow(10.,precdigits);
    if(roundr == 0) {
      std::cerr << "increase precision because edge weights too small" << std::endl;
      return -1;
    }
    printVal[idx]=-roundr;
    idx++;
    diag[u]+=roundr;
    diag[v]+=roundr;
  }

  printf("%%%%MatrixMarket matrix coordinate real symmetric\n%%\n%%Total Stretch %f\n", stretch);
  std::cout << n << ' ' << n  << ' ' << m+n-1+n << std::endl;

  
  for(int i=0; i < n-1; ++i) {
    std::cout << i+1 << ' ' << i+2 << ' ' << -1 << std::endl;
    
    diag[i]+=1;
    diag[i+1]+=1;
  }


  idx=0;
  for(iter=has.begin(); iter != has.end(); ++iter) {
    std::cout << (*iter).first+1 << ' ' << (*iter).second+1 << ' ' << printVal[idx] << std::endl;
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

