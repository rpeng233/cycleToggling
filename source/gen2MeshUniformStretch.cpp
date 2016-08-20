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

int a[5000][5000];
double ew[27000000];
std::map<double, int> mp;
double rS[27000000];

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
	
	int width=int(sqrt(n));
	n=width*width;
	
	int m = width*(width-1)*2-(n-1);

	int now=0;
	rep(i,0,width-1) a[0][i]=now, now++;
	rep(i,1,width-1) a[i][width-1]=now, now++;
	int flag=1;
	repd(i,width-1,1)
	{
		if (flag)
			repd(j,width-2,0)
				a[i][j]=now, now++;
		else
			rep(j,0,width-2)
				a[i][j]=now, now++;
		flag=1-flag;
	}
	
	
	cout << n << ' ' << m + n - 1 << ' ' << k << endl;

	rS[0]=0;
	rep(i,1,n) rS[i]=rS[i-1]+rand()%1000+1;
	
	for (size_t i = 0; i < n - 1; i++) {
		cout << i << ' ' << i + 1 << ' ' << rS[i+1]-rS[i] << endl;
	}
	
	
	size_t s = 0;
	rep(i,0,width-1)
		rep(j,0,width-2)
			if (abs(a[i][j]-a[i][j+1])!=1)
			{
				double w=abs(rS[a[i][j]]-rS[a[i][j+1]]); 
				cout << a[i][j] << ' ' << a[i][j+1] << ' ' << w << endl;
				ew[s]=w; s++;
			}
	
	rep(i,0,width-1)
		rep(j,0,width-2)
			if (abs(a[j][i]-a[j+1][i])!=1)
			{
				double w=abs(rS[a[j][i]]-rS[a[j+1][i]]); 
				cout << a[j][i] << ' ' << a[j+1][i] << ' ' << w << endl;
				ew[s]=w; s++;
			}
	assert(s==m);

	for (size_t i = 0; i < n - 1 + m; i++) {
		if (i < n - 1) {
			cout << "1.0" << endl;
		} else {
			cout << "0.0" << endl;
		}
	}
	double totw=0;
	rep(i,0,m-1) 
	{
		totw+=ew[i];
		mp[totw]=i;
	}
	
	/*for (size_t i = 0; i < k; ++i) {
		printf("%d\n", n + sample(totw,m));
		}*/

	return 0;
}
