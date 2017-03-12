#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>

int main(int argc, char *argv[]) {

  if(argc < 4) {
    std::cerr << "Please specify filename, n, and seed." << std::endl;
    return -1;
  }
  
  std::string rfile(argv[1]);
  std::string mmfile=rfile+".mtx";
  rfile=rfile+".txt";

  std::ofstream rfileout, mmfileout;
  rfileout.open(rfile);
  mmfileout.open(mmfile);

  size_t n = atol(argv[2]);
  size_t hop = 2;
 
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

  size_t seed = atol(argv[3]);

  srand(seed);
  int m=n-hop;
  
  mmfileout << "%%MatrixMarket matrix coordinate real symmetric" << std::endl;
  mmfileout << "%%" << std::endl;
  mmfileout << "%%Total Stretch " << m << std::endl;
  mmfileout << n << ' ' << n  << ' ' << m+n-1+n << std::endl;

  rfileout << n << ' ' << m+n-1 << std::endl;


  rS[0]=0;
  for(int i=0; i < n-1; ++i) {
    double r=rand()%1000+1;
    double roundr=round(pow(10.,precdigits)/r)/pow(10.,precdigits);
    if(roundr == 0 || isnan(roundr)) {
      std::cerr << "increase precision beyond " << precdigits
                << " because edge weights too small" << std::endl;
      return -1;
    }
    mmfileout << i+1 << ' ' << i+2 << ' ' << -roundr << std::endl; 
    rfileout << i << ' ' << i+1 << ' ' << r << std::endl; 
    diag[i]+=roundr;
    diag[i+1]+=roundr;
    rS[i+1]=rS[i]+r;
  }

  for(int i=0; i < n-hop; ++i) {
    int u=i;
    int v=i+hop;     
    double roundr=round(pow(10.,precdigits)/(rS[v]-rS[u]))/pow(10.,precdigits);
    if(roundr == 0 || isnan(roundr)) {
      std::cerr << "increase precision beyond " << precdigits
                << " because edge weights too small" << std::endl;
      return -1;
    }
    rfileout << u << ' ' << v << ' ' << rS[v]-rS[u] << std::endl;
    mmfileout << u+1 << ' ' << v+1 << ' ' << -roundr << std::endl;
    diag[u]+=roundr;
    diag[v]+=roundr;
  }



  for(int i=0; i < n; ++i) {
    mmfileout << i+1 << ' ' << i+1 << ' ' << std::setprecision(precdigits+1) << diag[i] << std::endl;
  }

  delete diag;
  delete rS;

  rfileout.close();
  mmfileout.close();

  return 0;
}

