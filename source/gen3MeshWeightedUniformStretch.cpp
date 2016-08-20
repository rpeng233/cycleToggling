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
	if (argc < 3) {
		cerr << "Please specify n and k." << endl;
		return -1;
	}

	size_t n = atol(argv[1]);	
	size_t k = atol(argv[2]);
	
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
	
	fprintf(stderr,"0 complete\n");
	
	cout << n << ' ' << m + n - 1 << ' ' << k << endl;

	for (size_t i = 0; i < n - 1; i++) {
		r[i+1] = rand()%1000+1;
		cout << i << ' ' << i + 1 << ' ' << r[i+1] << endl;
		r[i+1]+=r[i];
	}
	
	fprintf(stderr,"1 complete %d %d\n",n,m);
	
	size_t s = 0;
	rep(i,0,width-1)
		rep(j,0,width-1)
			rep(k,0,width-2)
			{
				if (abs(a[i][j][k]-a[i][j][k+1])!=1)
				{
					double w=abs(r[a[i][j][k]]-r[a[i][j][k+1]]); 
					cout << a[i][j][k] << ' ' << a[i][j][k+1] << ' ' << w << endl;
					ew[s]=w; s++;
				}
				if (abs(a[i][k][j]-a[i][k+1][j])!=1)
				{
					double w=abs(r[a[i][k][j]]-r[a[i][k+1][j]]); 
					cout << a[i][k][j] << ' ' << a[i][k+1][j] << ' ' << w << endl;
					ew[s]=w; s++;
				}
				if (abs(a[k][i][j]-a[k+1][i][j])!=1)
				{
					double w=abs(r[a[k][i][j]]-r[a[k+1][i][j]]); 
					cout << a[k][i][j] << ' ' << a[k+1][i][j] << ' ' << w << endl;
					ew[s]=w; s++; 
				}
			}
	
	assert(s==m);
	
	fprintf(stderr,"2 complete\n");
	

	for (size_t i = 0; i < n - 1 + m; i++) {
		if (i < n - 1) {
			cout << "1.0" << endl;
		} else {
			cout << "0.0" << endl;
		}
	}

	fprintf(stderr,"3 complete\n");
	
	double totw=0;
	rep(i,0,m-1) 
	{
		totw+=ew[i];
		mp[totw]=i;
	}
	/*
	for (size_t i = 0; i < k; ++i) {
		printf("%d\n", n + sample(totw,m));
		}*/

	return 0;
}
