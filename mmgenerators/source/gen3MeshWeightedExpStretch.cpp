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
  
  int width=int(pow(n,1.0/3.0));
  n=width*width*width;
  
  size_t precdigits=6;

  int*** a = new int**[width];
  for(int i=0; i < width; ++i) {
    a[i] = new int*[width];
    for(int j=0; j < width; ++j) {
      a[i][j] = new int[width];
    }
  }

  double *diag = new double[n];
  for(int i=0; i < n; ++i) {
    diag[i]=0;
  }

  int m=width*width*(width-1)*3-(n-1);
  
  double *rS = new double[n];

  double *printVal = new double[m];
  int idx=0;

  double stretch=0;

  int now=0;
  int flag=1, flag2=1;

  for(int k=0; k <= width-1; ++k) {
    if(flag2) {
      for(int i=0; i <= width-1; ++i) {
        if(flag) {
          for(int j=width-1; j >= 0; j--) {
            a[k][i][j]=now++;
          }
        }
        else {
          for(int j=0; j <= width-1; ++j) {
            a[k][i][j]=now++;
          }
        }
        flag=1-flag;
      }
    }
    else {
      for(int i=width-1; i >= 0; --i) {
        if(flag) {
          for(int j=width-1; j >= 0; --j) {
            a[k][i][j]=now++;
          }
        }
        else {
          for(int j=0; j <= width-1; ++j) {
            a[k][i][j]=now++;
          }
        }
        flag=1-flag;
      }
    }
    flag2=1-flag2;
  }

  rS[0]=0;
  for(int i=0; i < n-1; ++i) {
    rS[i+1]=rS[i]+rand()%1000+1;
  }


  for(int i=0; i <= width-1; ++i) {
    for(int j=0; j <= width-1; ++j) {
      for(int k=0; k <= width-2; ++k) {
        if(abs(a[i][j][k]-a[i][j][k+1])!=1) {
          double r=abs(rS[a[i][j][k]]-rS[a[i][j][k+1]]);
          while(double(rand())/double(RAND_MAX) < 0.5) {
            r/=2;
          }
          stretch+=abs(rS[a[i][j][k]]-rS[a[i][j][k+1]])/r;
          double roundr=round(pow(10.,precdigits)/r)/pow(10.,precdigits);
          printVal[idx]=-roundr;
          diag[a[i][j][k]]+=roundr;
          diag[a[i][j][k+1]]+=roundr;
          idx++;
        }
        if(abs(a[i][k][j]-a[i][k+1][j])!=1) {
          double r=abs(rS[a[i][k][j]]-rS[a[i][k+1][j]]); 
          while(double(rand())/double(RAND_MAX) < 0.5) {
            r/=2;
          }
          double roundr=round(pow(10.,precdigits)/r)/pow(10.,precdigits);
          stretch+=abs(rS[a[i][k][j]]-rS[a[i][k+1][j]])/r;
          printVal[idx]=-roundr;
          diag[a[i][k][j]]+=roundr;
          diag[a[i][k+1][j]]+=roundr;
          idx++;
        }
        if(abs(a[k][i][j]-a[k+1][i][j])!=1) {
          double r=abs(rS[a[k][i][j]]-rS[a[k+1][i][j]]);
          while(double(rand())/double(RAND_MAX) < 0.5) {
            r/=2;
          }
          double roundr=round(pow(10.,precdigits)/r)/pow(10.,precdigits);
          stretch+=(rS[a[k][i][j]]-rS[a[k+1][i][j]])/r;
          diag[a[k][i][j]]+=roundr;
          diag[a[k+1][i][j]]+=roundr;
          idx++;
        }
      }
    }
  }

  idx=0;
  printf("%%%%MatrixMarket matrix coordinate real symmetric\n%%\n%%Total Stretch %f\n", stretch);
  std::cout << n << ' ' << n  << ' ' << m+n-1+n << std::endl;


  for (int i=0; i < n-1; ++i) {
    double roundr = round(pow(10.,precdigits)/(rS[i+1]-rS[i]))/pow(10.,precdigits);
    std::cout << i+1 << ' ' << i+2 << ' ' << -roundr << std::endl;
    diag[i]+=roundr;
    diag[i+1]+=roundr;
  }
  
  for(int i=0; i <= width-1; ++i) {
    for(int j=0; j <= width-1; ++j) {
      for(int k=0; k <= width-2; ++k) {
        if(abs(a[i][j][k]-a[i][j][k+1])!=1) {
          std::cout << a[i][j][k]+1 << ' ' << a[i][j][k+1]+1 << ' ' << printVal[idx] << std::endl;
          idx++;
        }
        if(abs(a[i][k][j]-a[i][k+1][j])!=1) {
          std::cout << a[i][k][j]+1 << ' ' << a[i][k+1][j]+1 << ' ' << -printVal[idx] << std::endl;
          idx++;
        }
        if(abs(a[k][i][j]-a[k+1][i][j])!=1) {
          std::cout << a[k][i][j]+1 << ' ' << a[k+1][i][j]+1 << ' ' << -printVal[idx] << std::endl;
          idx++;
        }
      }
    }
  }


  for(int i=0; i < n; ++i) {
    std::cout << i+1 << ' ' << i+1 << ' ' << std::setprecision(precdigits+1) << diag[i] << std::endl;
  }

  for(int i=0; i < width; ++i) {
    for(int j=0; j < width; j++) {
      delete a[i][j];
    }
    delete a[i];
  }

  delete a;
  delete diag;
  delete rS;
  delete printVal;

  return 0;
}
