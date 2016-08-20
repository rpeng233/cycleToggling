#include <cstdlib>
#include <cstdio> 
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>
#include <map>

using std::cout;
using std::cerr;
using std::endl;

#define SIZE(x) (int((x).size()))
#define rep(i,l,r) for (int i=(l); i<=(r); i++)
#define repd(i,r,l) for (int i=(r); i>=(l); i--)
#define rept(i,c) for (typeof((c).begin()) i=(c).begin(); i!=(c).end(); i++)

int a[300][300][300];
double ew[27000000];
std::map<double, int> mp;
double r[27000000];
double diag[27000000];

double genrand()
{
	double r=RAND_MAX;
	return rand()/r/r+rand()/r;	//windows无人权
}

int sample(double totw, int m)
{
	double v=totw*genrand();
	std::map<double, int>::iterator it=mp.lower_bound(v);
	if (it==mp.end()) it--;
	return it->second;
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		cerr << "Please specify n." << endl;
		return -1;
	}

	size_t n = atol(argv[1]);	
		
	int width=int(pow(n,1.0/3.0));
	n=width*width*width;
	
	int m = width*width*(width-1)*3-(n-1);

	int now=0;
	int flag=1, flag2=1;
	rep(k,0,width-1)
	{
		if (flag2)
			rep(i,0,width-1)
			{
				if (flag)
					repd(j,width-1,0)
						a[k][i][j]=now, now++;
				else
					rep(j,0,width-1)
						a[k][i][j]=now, now++;
				flag=1-flag;
			}
		else
			repd(i,width-1,0)
			{
				if (flag)
					repd(j,width-1,0)
						a[k][i][j]=now, now++;
				else
					rep(j,0,width-1)
						a[k][i][j]=now, now++;
				flag=1-flag;
			}
		flag2=1-flag2;
	}

	printf("%%%%MatrixMarket matrix coordinate real symmetric\n%%\n");
	cout << n << ' ' << n  << ' ' << m+n+n << endl;


	for (size_t i = 0; i < n - 1; i++) {
		r[i+1] = 1;
		cout << i+1 << ' ' << i + 2 << ' ' << -1./r[i+1] << endl;
		diag[i]+=1./r[i+1];
		diag[i+1]+=1./r[i+1];
		r[i+1]+=r[i];
	}
	

	
	size_t s = 0;
	rep(i,0,width-1)
		rep(j,0,width-1)
			rep(k,0,width-2)
			{
				if (abs(a[i][j][k]-a[i][j][k+1])!=1)
				{
					double w=abs(r[a[i][j][k]]-r[a[i][j][k+1]]); 
					while (double(rand())/double(RAND_MAX)<0.5) w/=2;
					cout << a[i][j][k]+1 << ' ' << a[i][j][k+1]+1 << ' ' << -1./w << endl;
					diag[a[i][j][k]]+=1./w;
					diag[a[i][j][k+1]]+=1./w;
					ew[s]=w; s++;
				}
				if (abs(a[i][k][j]-a[i][k+1][j])!=1)
				{
					double w=abs(r[a[i][k][j]]-r[a[i][k+1][j]]); 
					while (double(rand())/double(RAND_MAX)<0.5) w/=2;
					cout << a[i][k][j]+1 << ' ' << a[i][k+1][j]+1 << ' ' << -1./w << endl;
					diag[a[i][j][j]]+=1./w;
					diag[a[i][k+1][j]]+=1./w;
					ew[s]=w; s++;
				}
				if (abs(a[k][i][j]-a[k+1][i][j])!=1)
				{
					double w=abs(r[a[k][i][j]]-r[a[k+1][i][j]]); 
					while (double(rand())/double(RAND_MAX)<0.5) w/=2;
					cout << a[k][i][j]+1 << ' ' << a[k+1][i][j]+1 << ' ' << -1./w << endl;
					diag[a[k][i][j]]+=1./w;
					diag[a[k+1][i][j]]+=1./w;
					ew[s]=w; s++; 
				}
			}
	rep(i,0,n-1)			
	  cout << i+1 << ' ' << i+1 << ' ' << diag[i] << std::endl;
	

	return 0;
}
