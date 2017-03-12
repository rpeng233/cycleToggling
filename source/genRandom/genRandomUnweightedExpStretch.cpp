#include <iostream>
#include <iomanip>
#include <math.h>
#include <set>
#include <fstream>

std::set<std::pair<int,int> > has;

int main(int argc, char *argv[]) {

  if(argc < 4) {
    std::cerr << "Please specify filename, n, and m" << std::endl;
  }
  
  std::string rfile(argv[1]);
  std::string mmfile=rfile+".mtx";
  rfile=rfile+".txt";

  std::ofstream rfileout, mmfileout;
  rfileout.open(rfile);
  mmfileout.open(mmfile);

  size_t n=atol(argv[2]);
  size_t m=atol(argv[3]);
  size_t precdigits=6;
  if(argc == 5) {
    precdigits=atol(argv[4]);
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
    printVal[idx]=newr;
    idx++;
  }

  mmfileout << "%%MatrixMarket matrix coordinate real symmetric" << std::endl;
  mmfileout << "%%" << std::endl;
  mmfileout << "%%Total Stretch " << stretch << std::endl;
  mmfileout << n << ' ' << n  << ' ' << m+n-1+n << std::endl;

  rfileout << n << ' ' << m+n-1 << std::endl;
  
  for(int i=0; i < n-1; ++i) {
    mmfileout << i+1 << ' ' << i+2 << ' ' << -1 << std::endl;
    rfileout << i << ' ' << i+1 << ' ' << 1 << std::endl;
    diag[i]+=1;
    diag[i+1]+=1;
  }


  idx=0;
  for(iter=has.begin(); iter != has.end(); ++iter) {
    double roundr=round(pow(10.,precdigits)/printVal[idx])/pow(10.,precdigits);
    if(roundr == 0) {
      std::cerr << "increase precision because edge weights too small" << std::endl;
      return -1;
    } 
    diag[(*iter).first]+=roundr;
    diag[(*iter).second]+=roundr;
    mmfileout << (*iter).first+1 << ' ' << (*iter).second+1 << ' ' << -roundr << std::endl;
    rfileout << (*iter).first << ' ' << (*iter).second << ' ' << printVal[idx] << std::endl;
    idx++;
  }
  
  for(int i=0; i < n; ++i) {
    mmfileout << i+1 << ' ' << i+1 << ' ' << std::setprecision(precdigits+1) << diag[i] << std::endl;
  }

  delete diag;
  delete rS;
  delete printVal;

  mmfileout.close();
  rfileout.close();

  return 0;
}

