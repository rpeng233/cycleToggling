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
#include "../util/io.h"
#include "../util/graph.h"

#define maxm 9000010	//edge #
#define maxo 10000010	//query #

static int n,m;
static bool init_run = true;

struct etype
{
	int which;
	double fPath, rPath;
};
static vector< pair<int,etype> > e[maxn];

namespace HLD
{
	struct atype
	{
		int parent, type;		//type=0: left child; 1: right child; 2: virtual edge
		double sum, sum_resistance, tag, rPath, fPath;
	};
	
	struct btype
	{
		int which, lc, rc;
		double rtag;
	};
	
	struct data_type
	{
		double rPath, fPath;
		int which, id;
	};
	
	static etype parente[maxn];
	static int parent[maxn], size[maxn], weight[maxn], isHeavy[maxn], isLeaf[maxn], s[maxn], tl[maxn];
	static btype d2[maxn];
	
	static int all;
	static data_type d[maxn];
	static atype ds[maxn];
	
	void dfs(int cur)
	{
		int maxs=-1, maxsi=0; size[cur]=1;
		rept(it,e[cur])
		{
			dfs(it->first); parent[it->first]=cur; parente[it->first]=it->second;
			size[cur]+=size[it->first];
			if (size[it->first]>maxs)
			{
				maxs=size[it->first]; maxsi=it->first;
			}
		}
		if (maxsi) isHeavy[maxsi]=1; else isLeaf[cur]=1;
		weight[cur]=size[cur];
		if (maxsi) weight[cur]-=size[maxsi];
	}
	
	int find_mid(int l, int r)
	{
		int target=(s[l-1]+s[r]+1)/2;
		while (l<=r)
		{
			int mid=(l+r)/2;
			if (s[mid-1]<target && target<=s[mid]) return mid;
			if (target<=s[mid-1]) r=mid-1; else l=mid+1;
		}
		assert(0);
	}
	
	pair<int,double> build(int l, int r)
	{
		int root=find_mid(l,r); double s=0;
		int id=d[root].id; ds[id].tag=0; ds[id].sum=0;
		all++; tl[all]=id;
		if (root>l) 
		{
			pair<int,double> lc=build(l,root-1); d2[id].lc=lc.first;
			ds[lc.first].parent=id; ds[lc.first].type=0; 
			ds[lc.first].rPath=d[root].rPath; ds[lc.first].fPath=0;
			d2[lc.first].which=d[root].which;
			ds[id].sum_resistance=d[root].rPath+lc.second;
		}
		else	
		{
			ds[id].sum_resistance=0;
			d2[id].lc=0;
		}
		if (root<r)
		{
			pair<int,double> rc=build(root+1,r); d2[id].rc=rc.first;
			ds[rc.first].parent=id; ds[rc.first].type=1; 
			ds[rc.first].rPath=d[root+1].rPath; ds[rc.first].fPath=0;
			d2[rc.first].which=d[root+1].which;
			s+=rc.second+d[root+1].rPath;
		}
		else	d2[id].rc=0;
		return make_pair(id,s+ds[id].sum_resistance);
	}

	void construct(int m, int head)
	{
		s[0]=0;
		rep(i,1,m) s[i]=s[i-1]+weight[d[i].id];
		pair<int,double> rt=build(1,m);
		ds[rt.first].parent=head; ds[rt.first].type=2; 
		ds[rt.first].rPath=d[1].rPath; ds[rt.first].fPath=0;
		d2[rt.first].which=d[1].which;
	}
	
	void hld()
	{
		memset(isHeavy,0,sizeof isHeavy);
		memset(isLeaf,0,sizeof isLeaf);
		dfs(1);
		all=0;
		rep(i,1,n)
			if (isLeaf[i])
			{
				int ci=0, now=i;
				while (isHeavy[now])
				{
					ci++; now=parent[now];
				}
				ci++; int z=parent[now]; now=i;
				repd(k,ci,1)
				{
					d[k].rPath=parente[now].rPath; 
					d[k].fPath=parente[now].fPath; 
					d[k].which=now;
					d[k].id=now; 
					now=parent[now];
				}
				construct(ci,z);
			}
	}
	
	double query(int x)
	{
		double ans=0;
		while (x)
		{
			int flag=1; double rp=0;
			while (1)
			{
				if (flag) 
				{
					ans+=ds[x].sum; 
					rp+=ds[x].sum_resistance;
				}
				else	
				{
					ans+=rp*ds[x].tag;
				}
				if (ds[x].type!=0) 
				{
					ans+=ds[x].rPath*ds[x].fPath;
					rp+=ds[x].rPath;
				}
				if (ds[x].type==0) flag=0; else flag=1;
				if (ds[x].type==2) break;
				x=ds[x].parent;
			}
			x=ds[x].parent;
		}
		return ans;
	}
	
	void serere(int x, double value)
	{
		while (x)
		{
			int flag=1; double rp=0;
			while (1)
			{
				if (flag) 
				{
					ds[x].tag+=value;
					double tmp=ds[x].sum_resistance*value;
					ds[x].sum+=tmp; rp+=tmp;
				}
				else	
				{
					ds[x].sum+=rp;
				}
				if (ds[x].type!=0) 
				{
					ds[x].fPath+=value;
					rp+=ds[x].rPath*value;
				}
				if (ds[x].type==0) flag=0; else flag=1;
				if (ds[x].type==2) break;
				x=ds[x].parent;
			}
			x=ds[x].parent;
		}
	}
	
	void dump(double *ans)
	{
		rep(i,1,n) ans[i]=0; 
		rep(i,1,n) d2[i].rtag=0;
		rep(tx,1,all)
		{
			int i=tl[tx];
			if (ds[i].type!=0) ans[d2[i].which]+=ds[i].fPath;
			if (d2[i].lc!=0) 
			{
				ans[d2[d2[i].lc].which]+=ds[i].tag+d2[i].rtag;
				d2[d2[i].lc].rtag+=ds[i].tag+d2[i].rtag;
			}
			if (d2[i].rc!=0)
			{
				d2[d2[i].rc].rtag+=d2[i].rtag;
				ds[d2[i].rc].fPath+=d2[i].rtag;
			}
		}
	}
}

namespace LCA
{
	static int depth[maxn]; 
	static int p[maxn][22];
	static int lg2[maxn];
	
	void dfs(int cur, int pre, int dep)
	{
		depth[cur]=dep; p[cur][0]=pre;
		rep(i,1,lg2[dep]) p[cur][i]=p[p[cur][i-1]][i-1];
		rept(it,e[cur]) dfs(it->first,cur,dep+1);
	}
		
	void init()
	{
		lg2[1]=0;
		rep(i,2,maxn-1) lg2[i]=lg2[i>>1]+1;
		dfs(1,0,0);
	}
	
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
}

static double resistance[maxm], f[maxm], b[maxn], tfPath[maxn];
static int ote[maxm][2];
static double energy;

static int kParam;
static int tt;
static double bnorm;

double sqr(double x) { return x*x; }

namespace Solver
{
	static int all;
	static double osum_r[maxn], osum_rf[maxn];
	static int lca[maxm], lis[maxn];
	
	void dfs(int cur)
	{
		all++; lis[all]=cur;
		rept(it,e[cur])
		{
			osum_r[it->first]=osum_r[cur]+it->second.rPath;
			osum_rf[it->first]=osum_rf[cur]+it->second.rPath*it->second.fPath;
			dfs(it->first);
		}
	}
	
	void init()
	{
		HLD::hld();
		LCA::init();
		osum_r[1]=0; osum_rf[1]=0; all=0; dfs(1);
		rep(i,1,m) lca[i]=LCA::lca(ote[i][0],ote[i][1]);
	}
	
	void toggle(int id)
	{
		int x=ote[id][0], y=ote[id][1], z=lca[id];
		double sum_rf=f[id]*resistance[id];
		double sum_r=resistance[id];
		sum_rf+=osum_rf[x]-osum_rf[y];
		sum_r+=osum_r[x]+osum_r[y]-osum_r[z]*2;
		sum_rf+=HLD::query(x)-HLD::query(y);
		double delta=sum_rf/sum_r;
		//mexPrintf("Toggle %d %d %d\n%.3lf %.3lf %.3lf\n",id,x,y,sum_rf,sum_r,delta);
		f[id]-=delta;
		HLD::serere(x,-delta); HLD::serere(y,delta);
		energy -= delta * sum_rf;
	}
	
	void dump(double *ans)
	{
		HLD::dump(ans);
	}
	
	void reset()
	{
		osum_rf[1]=0;
		rep(kk,1,n)
		{
			int i=lis[kk];
			if (i>1) osum_rf[i]=osum_rf[HLD::parent[i]]+tfPath[i]*HLD::parente[i].rPath;
		}
		
		rep(i,1,n)
		{
			HLD::ds[i].sum=0; HLD::ds[i].tag=0;
			HLD::ds[i].fPath=0;
		}
	}
	
	static double x[maxn], b1[maxn];
	
	int calculate_residual()
	{
		tt++;
		double eFlow = 0;
		double eVoltage = 0;
		rep(i,1,n) b1[i]=0;
		x[1] = 0;
		rep(kk,1,n)
		{
			int i=lis[kk];
			if (i==1) continue;
			double vv=tfPath[i]*HLD::parente[i].rPath;
			x[i] = x[HLD::parent[i]] + vv;
			eFlow += HLD::parente[i].rPath * sqr(tfPath[i]);
			eVoltage += sqr(vv) / HLD::parente[i].rPath;
			b1[HLD::parent[i]]-=tfPath[i];
			b1[i]+=tfPath[i];
		}
		rep(i,1,m)
		{
			eFlow += resistance[i] * sqr(f[i]);
			eVoltage += sqr(x[ote[i][0]] - x[ote[i][1]]) / resistance[i];
			double fx=(x[ote[i][0]] - x[ote[i][1]]) / resistance[i];
			b1[ote[i][0]]+=fx;
			b1[ote[i][1]]-=fx;
		}
		double xb = 0;
		rep(i,0,n-1) xb += x[i+1] * b[i];

		double gap = eFlow - (2 * xb - eVoltage);
		fprintf(stderr, "iteration %d: energy = %lf = %lf, xb = %lf, voltage energy = %lf, gap = %lf\n", tt, energy, eFlow, xb, eVoltage, gap);
		
		double b1norm = 0;
		rep(i,0,n-1) b1norm+=sqr(b1[i+1]-b[i]);
		b1norm=sqrt(b1norm);
                //printf("residual %0.6e\n",b1norm/bnorm);
		return (b1norm <= double(1e-5) * bnorm) || ( LL(tt) > 100000);
	}
	
	double calc_stretch(int id)
	{
		int u=ote[id][0], v=ote[id][1], w=lca[id];
		double rs=osum_r[u]+osum_r[v]-osum_r[w]*2;
		return rs/double(resistance[id]);
	}
}

static int order[maxo], tu[maxn], tv[maxn];
static double ans[maxn], trPath[maxn];

static double scale;
static int order_num;
static double sf[maxn], sw[maxm];

int sampleToggle()
{
	double x=double(rand())/double(RAND_MAX)*sw[m];
	int t=lower_bound(sw+1,sw+m+1,x)-sw;
	return t;
}

int main(int argc, char *argv[])
{
	if (argc!=3) { printf("Usage: ./linesolver graphFile rhsFile\n"); return 0; }
	
	string graphFile=argv[1];
	string rhsFile=argv[2];
	
	GraphSP g=IO::readGraphSP(graphFile);
	n=g.n;
	
	int all=0;
	rep(i,1,n)
		rept(it,g.e[i])
		{
			all++;
			tu[all]=i; tv[all]=it->first; trPath[all]=it->second;
		}
		
	m=0;
	rept(it,g.o)
	{
		m++; 
		ote[m][0]=get<0>(*it); ote[m][1]=get<1>(*it); 
		resistance[m]=get<2>(*it);
	}
	
	rep(i,1,n-1) tfPath[i]=0;
	rep(i,1,m) f[i]=0;
	
	rep(i,1,n-1)
	{
		etype t; t.rPath=trPath[i]; t.fPath=tfPath[i]; t.which=i;
		//if (tu[i]>tv[i]) swap(tu[i],tv[i]);
		e[tu[i]].push_back(make_pair(tv[i],t));
	}
	
	Solver::init();
		
	sw[0]=0; rep(i,1,m) sw[i]=sw[i-1]+(1+Solver::calc_stretch(i));
	
	order_num = n*4;
	
	Vec ds=IO::readMMVec(rhsFile);
	assert(ds.n==n);
	
	rep(i,0,n-1) b[i]=ds[i];
	
	double s = 0;
	rep(i,0,n-1)
	{
		s += b[i];
	}
	s = s / double(n);
	rep(i,0,n-1) b[i] -= s;
	
	bnorm = 0;
	rep(i,0,n-1) bnorm+=sqr(b[i]);
	bnorm = sqrt(bnorm);

	energy = 0;
	rep(i,1,n) sf[i]=0;
	repd(kk,n,1)
	{
		int i=Solver::lis[kk];
		sf[HLD::parent[i]]+=sf[i]-b[i-1];
		tfPath[i] = -(sf[i]-b[i-1]);
		energy += HLD::parente[i].rPath * sqr(tfPath[i]);
	}
		
	rep(i,1,m) f[i]=0;
	
	//sprintf(msg, ":::initial energy = %lf\n", energy);
	//mexWarnMsgIdAndTxt("MyToolbox:recursiveLineMex:lol", msg);
	
	static int tmp[maxo];
	
	double tsamplecost=0;
	clock_t t_start = clock();
	while (1)
	{
		Solver::reset();
		
		clock_t t_start2 = clock();
		rep(i,1,order_num) tmp[i]=sampleToggle();
		clock_t t_end2 = clock();
		tsamplecost+=double(t_end2 - t_start2) / double(CLOCKS_PER_SEC);
		
		rep(i,1,order_num) Solver::toggle(tmp[i]);

		Solver::dump(ans);
		
		rep(i,2,n) tfPath[i]+=ans[i];
		
		if (Solver::calculate_residual()) break;
	}
	clock_t t_end = clock();
	double tcost=double(t_end - t_start) / double(CLOCKS_PER_SEC);
	
	Vec vs(n);
	rep(i,0,n-1) vs[i]=Solver::x[i+1];
	Mat A=IO::constructMatrixFromGraph(g);
	double relres=(A*vs-ds).norm()/ds.norm();
	//printf("%.16lf\n",relres);
	//assert(relres<(1e-5)+(1e-9));
	printf("treesolver %0.6e %lld %.3lf %.3lf\n",relres,LL(tt)*order_num, tcost, tsamplecost);
	
	return 0;
}
