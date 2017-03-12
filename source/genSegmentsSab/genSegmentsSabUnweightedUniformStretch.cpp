#include <iostream>
#include <random>
#include <iomanip>
#include <math.h>
#include <set>
#include <fstream>

std::set<std::pair<int,int> > has;

int main(int argc, char *argv[])
{
  if (argc < 6) {
    std::cerr << "Please specify filename, n, m, seg, and sab resistance." << std::endl;
    return -1;
  }

  std::string rfile(argv[1]);
  std::string mmfile=rfile+".mtx";
  rfile=rfile+".txt";

  std::ofstream rfileout, mmfileout;
  rfileout.open(rfile);
  mmfileout.open(mmfile);

  size_t n=atol(argv[2]);
  size_t m=atol(argv[3]);
  size_t seg=atol(argv[4]);
  size_t sab_res=atol(argv[5]);
  
  size_t precdigits=6;
  
  int actualm = 0;
  int idx=0;
  size_t seg_size = n/seg;
  double *diag = new double[n];
  for(int i=0; i < n; ++i) {
    diag[i]=0;
  }

  double *rS = new double[n];
  double roundsab=round(pow(10.,precdigits)/sab_res)/pow(10.,precdigits);

  mmfileout << "%%MatrixMarket matrix coordinate real symmetric" << std::endl;
  mmfileout << "%%" << std::endl;
  mmfileout << "%%Total Stretch " << m << std::endl;
  mmfileout << n << ' ' << n  << ' ' << m+n-1+n << std::endl;

  rfileout << n << ' ' << m+n-1 << std::endl;


  rS[0]=0;
    
  for(int i=0; i < n-1; ++i) {
    if((i+1)%seg_size == 0 || (i+2)%seg_size == 0) {
      diag[i]+=roundsab;
      diag[i+1]+=roundsab;
      rS[i+1]=rS[i]+sab_res;
      mmfileout << i+1 << ' ' << i+2 << ' ' << -roundsab << std::endl;
      rfileout << i << ' ' << i+1 << ' ' << sab_res << std::endl;
    }
    else {
      diag[i]+=1;
      diag[i+1]+=1;
      rS[i+1]=rS[i]+1;
      mmfileout << i+1 << ' ' << i+2 << ' ' << -1 << std::endl;
      rfileout << i << ' ' << i+1 << ' ' << 1 << std::endl;
    }
  }

  

  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_int_distribution<int> uni(0,seg-1);
  std::uniform_int_distribution<int> uni2(0,seg_size-2);
  std::uniform_int_distribution<int> uni3(0,n-seg*seg_size-1);

  int segment;

  has.clear();

  for(int i=0; i < m; ++i) {
    int u=-1;
    int v=-1;
    
    while(u == v || v == u+1 || has.find(std::make_pair(u,v)) != has.end()) {

      //select segment
      segment=uni(rng);
  
      if(segment==seg-1 && (n-seg*seg_size) > 2) {
        u=uni3(rng)+segment*seg_size;
        v=uni3(rng)+segment*seg_size;
      }
      else {
        u=uni2(rng)+segment*seg_size;
        v=uni2(rng)+segment*seg_size;
      }

      if(u > v) {
        std::swap(u,v);
      }
    }
    
    has.insert(std::make_pair(u,v));

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
