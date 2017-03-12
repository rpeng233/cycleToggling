#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <cstring>
#include "../util/io.h"
#include "../util/graph.h"

using namespace std;


#define CHUNK 4
#define maxm 3000010
#define maxq 1000010

struct atype
{
	int p, x, l, r;		//parent, self and range of childs
	int *q;			//list of queries
	double deltaf, deltarf;	//tmp delta flow and voltage tag from this node to root
	int flag;			//tmp flag for contraction
	double resist, rfsum;		//SUM of voltage and resist from this node to root 
};

struct otetype
{
	int x,y;
	double resist, sum_r, flow;
};

struct etype
{
	int id;
	double resist, flow;
	etype() {}
	etype(int id, double resist, double flow): id(id), resist(resist), flow(flow) {}
};

namespace mempool
{
	atype a[4*maxq+maxn+1000];
	int ptr=0;
	
	atype *allocate(int sz)	//index 1~sz
	{
		int optr=ptr; ptr+=sz;
		return a+optr-1;
	}
	
	void restore(atype *pt)	//free everything allocated after pt (including pt itself) 
	{
		ptr=(pt-a)+1;
	}
}

int order[maxq], where[maxn];
otetype ote[maxm];

namespace Solver
{
	pair<atype*, int> contract(atype *a, int n, int ql, int qr)
	{
		atype *b=mempool::allocate(n);
		repd(i,n,1)
		{
			//flag: 0 not exist, >0 critical child id. real node has flag==i
			int flag=0;
			if (abs(*(a[i].q))<=qr) flag=i;
			if (flag==0)
			{
				int cnt=0;
				rep(j,a[i].l,a[i].r)
					if (a[j].flag>0)
					{
						cnt++; flag=a[j].flag;
					}
				if (cnt>=2 || i==1) flag=i;
			}
			a[i].flag=flag;
		}
		static int q[maxn][2];
		int head=1, tail=2, all=1; q[head][0]=1; q[head][1]=1;
		b[all].p=-1; b[all].q=a[1].q; b[all].rfsum=a[1].rfsum; b[all].resist=a[1].resist;
		while (head<tail)
		{
			int x=q[head][0], y=q[head][1]; head++;
			b[x].l=all+1; b[x].x=y;
			rep(i,a[y].l,a[y].r)
				if (a[i].flag>0)
				{
					int j=a[i].flag;
					all++; b[all].p=x; b[all].q=a[j].q; b[all].rfsum=a[j].rfsum; b[all].resist=a[j].resist; b[all].deltaf=0;
					q[tail][0]=all; q[tail][1]=j; tail++;
				}
			b[x].r=all;
		}
		mempool::restore(b+all);
		return make_pair(b,all);
	}

	void propogate(atype *a, int n, atype *b, int m, int calculate_new_voltage)
	{
		rep(i,1,n) { a[i].deltarf=a[i].deltaf; a[i].deltaf=0; }	//here deltarf is used as a tmp variable..
		rep(i,1,m) 
		{
			int x=b[i].x, y=b[i].p;
			a[x].q=b[i].q; a[x].deltaf+=b[i].deltaf; 
			if (y>0) a[b[y].x].deltaf-=b[i].deltaf;
		}
		repd(i,n,2) a[a[i].p].deltaf+=a[i].deltaf;
		if (calculate_new_voltage)
		{
			a[1].deltarf=0;
			rep(i,2,n) 
			{
				double t=a[a[i].p].deltarf+a[i].deltaf*(a[i].resist-a[a[i].p].resist);
				a[i].deltaf+=a[i].deltarf; a[i].deltarf=t; a[i].rfsum+=t; 
			}
		}
		else
		{
			rep(i,2,n) a[i].deltaf+=a[i].deltarf;
		}
	}

	double df[maxn];
	int tq[maxq][2];
	
	void brute_force(atype *a, int n, int ql, int qr)
	{
		rep(i,1,n) 
			while (abs(*(a[i].q))<=qr) 
			{
				int x=(*a[i].q);
				a[i].q++;
				if (x<0) tq[-x-ql][0]=i; else tq[x-ql][1]=i;
			}
		rep(i,1,n) a[i].deltaf=0;
		rep(i,ql,qr) 
		{
			int id=order[i];
			int x=tq[i-ql][0], y=tq[i-ql][1]; double sum_r=ote[id].sum_r;
			double sum_rf=ote[id].flow*ote[id].resist;
			sum_rf+=a[x].rfsum-a[y].rfsum;
			int t=x; while (t!=1) { sum_rf+=a[t].deltaf*(a[t].resist-a[a[t].p].resist); t=a[t].p; }
			t=y; while (t!=1) { sum_rf-=a[t].deltaf*(a[t].resist-a[a[t].p].resist); t=a[t].p; }
			double delta=sum_rf/sum_r;
			//printf("Toggle %d: %.3lf %.3lf %.3lf\n",id,sum_rf,sum_r,delta);
			ote[id].flow-=delta; df[ote[id].x]-=delta; df[ote[id].y]+=delta;
			while (x!=1) { a[x].deltaf-=delta; x=a[x].p; }
			while (y!=1) { a[y].deltaf+=delta; y=a[y].p; }
		}
	}

	void dc(atype *a, int n, int ql, int qr)
	{
		//printf("Tree @ %d %d:\n",ql,qr);
		//rep(i,1,n) printf("%d %d %d %.3lf %.3lf\n",a[i].p,i,a[i].x,a[i].resist,a[i].rfsum);
		//returns correct deltaf in a
		//nodes in a are stored in bfs order
		if (n<=100 || qr-ql<=100) { brute_force(a, n, ql, qr); return; }
		int len=(qr-ql+1)/CHUNK;
		rep(i,1,CHUNK)
		{
			int qm=ql+len-1; if (i==CHUNK) qm=qr;
			pair<atype*, int> z=contract(a, n, ql, qm);
			dc(z.first, z.second, ql, qm);
			propogate(a, n, z.first, z.second, (i<CHUNK));
			mempool::restore(z.first);
			ql=qm+1;
		}
	}
	
	double tmp[maxn];
	
	double *solve(atype *a, int n, int qn)
	{
		memset(df,0,sizeof df);
		dc(a,n,1,qn);
		rep(i,1,n) tmp[where[i]]=df[i];
		repd(i,n,2) tmp[a[i].p]+=tmp[i];
		rep(i,1,n) df[a[i].x]=tmp[i];
		return df;
	}
}

int p[maxn][30], depth[maxn], lg2[maxn];

int lca(int x, int y)
{
	if (depth[x]>depth[y]) swap(x,y);
	int z=depth[y]-depth[x];
	while (z)
	{
		y=p[y][lg2[z&-z]]; z-=z&-z;
	}
	repd(i,lg2[depth[x]],0)
		if (p[x][i]!=p[y][i])
		{
			x=p[x][i]; y=p[y][i];
		}
	if (x!=y) return p[x][0]; else return x;
}

vector<etype> e[maxn];

double sw[maxm], b[maxn], b1[maxn], x[maxn], sf[maxn];
int n,m;

int sampleToggle()
{
	double x=double(rand())/double(RAND_MAX)*sw[m];
	int t=lower_bound(sw+1,sw+m+1,x)-sw;
	return t;
}

void lemon(string filename, string gradname)
{
	//puts("Reading input");
	static int tu[maxn], tv[maxn], wtv[maxn];
	static double trPath[maxn], tfPath[maxn];
	FILE* ff = fopen(filename.c_str(), "r");
	int order_num; fscanf(ff,"%d%d", &n, &m); order_num=n;
	char buff[100];
	fgets(buff,100,ff);	//skip the "order_num" parameter
	rep(i,1,n-1)
	{
		fscanf(ff,"%d%d%lf", &tu[i], &tv[i], &trPath[i]);   
		tu[i]++; tv[i]++; wtv[tv[i]]=i;
	}
	m -= n - 1;
	rep(i,1,m) 
	{
		fscanf(ff,"%d%d%lf", &ote[i].x, &ote[i].y, &ote[i].resist);
		ote[i].x++; ote[i].y++;
	}
	
	rep(i,1,n-1) e[tu[i]].push_back(etype(tv[i],trPath[i],0));
	
	lg2[1]=0; rep(i,2,maxn-1) lg2[i]=lg2[i>>1]+1;
	static int q[maxn];
	int head=1, tail=2; q[1]=1; depth[1]=0; 
	atype *a=mempool::allocate(n);
	a[1].p=0; a[1].x=1; a[1].resist=0; a[1].rfsum=0; where[1]=1;
	while (head<tail)
	{
		int cur=q[head]; int dep=depth[cur];
		rep(i,1,lg2[dep]) p[cur][i]=p[p[cur][i-1]][i-1];
		a[head].l=tail;
		rept(it,e[cur])
		{
			depth[it->id]=dep+1; q[tail]=it->id; where[it->id]=tail; p[it->id][0]=cur;
			a[tail].p=head; a[tail].x=it->id; a[tail].deltaf=0;
			a[tail].resist=a[head].resist+it->resist; a[tail].rfsum=a[head].rfsum+it->resist*it->flow;
			tail++;
		}
		a[head].r=tail-1;
		head++;
	}
	
	sw[0]=0;
	rep(i,1,m)
	{
		int x=ote[i].x, y=ote[i].y; int z=lca(x,y);
		double pathrsum=a[where[x]].resist+a[where[y]].resist-2*a[where[z]].resist;
		ote[i].sum_r=ote[i].resist+a[where[x]].resist+a[where[y]].resist-2*a[where[z]].resist;
		sw[i]=sw[i-1]+pathrsum/ote[i].resist+1;
	}
	
	Vec ds=IO::readMMVec(gradname);
	assert(ds.n==n);
	
	rep(i,0,n-1) b[i]=ds[i];
	
	double bnorm = 0;
	rep(i,0,n-1) bnorm+=b[i]*b[i];
	bnorm = sqrt(bnorm);
	
	rep(i,1,n) sf[i]=0;
	repd(kk,tail-1,1)
	{
		int i=a[kk].x;
		int p=a[a[kk].p].x;
		sf[p]+=sf[i]-b[i-1];
		tfPath[wtv[i]] = -(sf[i]-b[i-1]);
	}
	
	rep(i,1,m) ote[i].flow=0;
	
	clock_t t_start = clock();
	double tsamplecost=0; int tt=0;
	while (1)
	{
		tt++;
		clock_t t_start2 = clock();
		static vector<int> qlist[maxn];
		rep(i,1,n) qlist[i].clear();
		rep(i,1,order_num) 
		{
			order[i]=sampleToggle();
			qlist[ote[order[i]].x].push_back(-i);
			qlist[ote[order[i]].y].push_back(i);
		}
		rep(i,1,n) qlist[i].push_back(order_num+1);
		int all=0;
		static int *whereq[maxn], qs[maxq*2+maxn];
		rep(i,1,n)
		{
			whereq[i]=qs+all;
			rept(it,qlist[i]) qs[all]=(*it), all++;
		}
		rep(i,1,tail-1) a[i].q=whereq[a[i].x];
		
		rep(i,2,tail-1) a[i].rfsum=a[a[i].p].rfsum+trPath[wtv[a[i].x]]*tfPath[wtv[a[i].x]];
		clock_t t_end2 = clock();
		tsamplecost+=double(t_end2 - t_start2) / double(CLOCKS_PER_SEC);
		
		double *ans=Solver::solve(a,n,order_num);
		rep(i,2,n) tfPath[wtv[i]]+=ans[i];
		
		x[1]=0;
		rep(i,1,n) b1[i]=0;
		rep(i,2,tail-1)
		{
			x[a[i].x]=x[a[a[i].p].x]+trPath[wtv[a[i].x]]*tfPath[wtv[a[i].x]];
			b1[a[a[i].p].x]-=tfPath[wtv[a[i].x]];
			b1[a[i].x]+=tfPath[wtv[a[i].x]];
		}
		rep(i,1,m)
		{
			double fx=(x[ote[i].x] - x[ote[i].y]) / ote[i].resist;
			b1[ote[i].x]+=fx;
			b1[ote[i].y]-=fx;
		}
		
		double b1norm = 0;
		rep(i,0,n-1) b1norm+=(b1[i+1]-b[i])*(b1[i+1]-b[i]);
		b1norm=sqrt(b1norm);
		
		double relres=b1norm/bnorm;
                fprintf(stderr, "residual %0.6e\n",relres);
		if ((relres<1e-5) || (LL(tt) > 100000)) break;
	}
	
	clock_t t_end = clock();
	double tcost=double(t_end - t_start) / double(CLOCKS_PER_SEC);
	
	GraphSP g=IO::readGraphSP(filename);
	
	Vec vs(n);
	rep(i,0,n-1) vs[i]=x[i+1];
	Mat A=IO::constructMatrixFromGraph(g);
	double relres=(A*vs-ds).norm()/ds.norm();
	//printf("%.16lf\n",relres);
	//assert(relres<(1e-5)+(1e-9));
	printf("rectreesolver %0.6e %lld %.3lf %.3lf\n",relres,LL(tt)*order_num, tcost, tsamplecost);
}

int main(int argc, char *argv[])
{
	if (argc!=3) { printf("Usage: ./linesolver graphFile rhsFile\n"); return 0; }
	ios::sync_with_stdio(true);
	lemon(argv[1],argv[2]);
	return 0;
}

