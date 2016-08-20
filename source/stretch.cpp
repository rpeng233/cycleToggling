#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <cstring>
#include <cassert>
#include "io.h"
#include "graph.h"

/**
 * please adjust bitn so that 2^bitn > (n+10)/block (do NOT ignore +10)
 * maxl > m+10
 */

#define block 20

int bitn, offset;

#define maxm 3000010
#define maxl 4000010

bool init_run = true;


struct data_type
{
	double fPath, rPath, ts;
};
data_type d[maxl];

double resistance[maxl];
double rS[maxn], sw[maxm];
int e[maxl][2];

int main(int argc, char *argv[])
{
	if (argc!=2) { printf("Usage: ./stretch graphFile"); return 0; }
	
	string graphFile=argv[1];
	
	GraphSP g=IO::readGraphSP(graphFile);
	int n=g.n;
	
	int z=1; bitn=0;
	while (z<=n/block+10) 
	{
		z*=2; bitn++;
	}
	offset = 1<<bitn;
		
	rS[1]=0;
	rep(i,1,n-1) 
	{
		assert(g.e[i].size()==1);
		assert(g.e[i][0].first==i+1);
		d[i].rPath=g.e[i][0].second;
		rS[i+1]=rS[i]+g.e[i][0].second;
	}
	
	int m=0;
	sw[0]=0;
	rept(it,g.o)
	{
		m++; 
		e[m][0]=get<0>(*it)-1; e[m][1]=get<1>(*it)-1; 
		if (e[m][0]>e[m][1]) swap(e[m][0],e[m][1]);
		resistance[m]=get<2>(*it);
		sw[m]=sw[m-1]+(rS[e[m][1]+1]-rS[e[m][0]+1])/resistance[m]+1;
	}
	printf("%.16lf\n",sw[m]);
	
	
	return 0;
}
