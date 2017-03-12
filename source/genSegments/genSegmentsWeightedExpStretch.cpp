#include <iostream>
#include <random>
#include <iomanip>
#include <math.h>
#include <set>
#include <fstream>

std::set<std::pair<int,int> > has;

int main(int argc, char *argv[])
{
  if (argc < 5) {
    std::cerr << "Please specify filename, n, m, and seg." << std::endl;
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
  
  size_t precdigits=6;

  if(argc==6) {
    precdigits=atol(argv[5]);
  }
  
  size_t seg_size = n/seg;
  double *diag = new double[n];
  for(int i=0; i < n; ++i) {
    diag[i]=0;
  }

  double *rS = new double[n];

  double *printVal = new double[m];

  rS[0]=0;
  for(int i=0; i < n-1; ++i) {
    rS[i+1]=rS[i]+rand()%1000+1;
  }
  
  double stretch=0;
  int idx=0;

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
  }

  std::set<std::pair<int,int> >::iterator iter;
  for(iter=has.begin(); iter != has.end(); ++iter) {
    int u=(*iter).first;
    int v=(*iter).second;
    
    double r=rS[v]-rS[u];
    double newr=r;
    double maxr=0;

    for(int l=0; l < v-u; ++l) {
      double tempr=rS[u+l+1]-rS[u+l];
      if(tempr > maxr) {
        maxr=tempr;
      }
    }
    
    while (double(rand())/double(RAND_MAX) < 0.5) {
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
    printVal[idx]=newr;
    idx++;
  }

  mmfileout << "%%MatrixMarket matrix coordinate real symmetric" << std::endl;
  mmfileout << "%%" << std::endl;
  mmfileout << "%%Total Stretch " << stretch << std::endl;
  mmfileout << n << ' ' << n  << ' ' << m+n-1+n << std::endl;

  rfileout << n << ' ' << m+n-1 << std::endl;

  for(int i=0; i < n-1; ++i) {
    double roundr=round(pow(10.,precdigits)/(rS[i+1]-rS[i]))/pow(10.,precdigits);
    if(roundr == 0 || isnan(roundr)) {
      std::cerr << "increase precision beyond " << precdigits
                << " because edge weights too small" << std::endl;
      return -1;
    }
    mmfileout << i+1 << ' ' << i+2 << ' ' << -roundr << std::endl;
    rfileout << i << ' ' << i+1 << ' ' << rS[i+1]-rS[i] << std::endl;
    diag[i]+=roundr;
    diag[i+1]+=roundr;
  }

  idx=0;
  for(iter=has.begin(); iter != has.end(); ++iter) {
    double roundr=round(pow(10.,precdigits)/printVal[idx])/pow(10.,precdigits);
    if(roundr == 0 || isnan(roundr)) {
      std::cerr << "increase precision beyond " << precdigits
                << " because edge weights too small" << std::endl;
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

  rfileout.close();
  mmfileout.close();

  return 0;
}
