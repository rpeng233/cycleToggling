/********************************************************************************
 * 
 * KOSZ Algorithm
 * 
 * Warning: the graph is 1-INDEXED, all inputs are RESISTANCE, not WEIGHT
 * 
 * Public API:
 * 	KOSZ(g): Constructor. 
 * 		GraphSP g: the graph 
 * 
 * 	(AbstractSolver):	implicit typecast to AbstractSolver. The initial guess x0 will be ignored.
 * 
 * 	setb(Vec b):
 * 		Sets the RHS of the equation to b
 * 
 * 	int sampleOTE():
 * 		Samples an off-tree edge w.p. propotional to 1+stretch
 * 		Returns an int k in [0,m), where m is total # of off-tree edges
 * 		The sampled off-tree edge is g.o[k]
 * 
 * 	solver:
 * 		... = solve()
 * 		... = solve(tol)
 * 		... = solve(tol, maxit)
 * 		... = solve(tol, maxit, batchnum)
 * 		LHS accepts any formats defined in abstractsolver.h
 * 
 * 	Inputs:
 * 		tol:		desired tolerance (default 1e-6)
 * 		maxit:	max # of batches to perform, default -1 (unlimited)
 * 		batchnum:	# of toggles in a batch, default -1 (4*n)
 * 	
 * 	Outputs:
 * 		see abstractsolver.h
 * 
 * Examples:
 * 	GraphSP g; Vec b;
 * 	KOSZ s1=KOSZ(g);			//construct as KOSZ Solver
 * 	Vec x=s1.solve(b);
 * 	AbstractSolver s2=KOSZ(g);	//construct as AbstractSolver
 * 	Vec x=s2.solve(b);
 * 
 ********************************************************************************/

#ifndef __KOSZ_H__
#define __KOSZ_H__

#include "common.h"
#include "matrix.h"
#include "abstractsolver.h"
#include "treefinder.h"

struct KOSZ
{
	int n,m;
	
	struct DS
	{
		int n;
		struct nodetype
		{
			int parent, type;		//type=0: left child; 1: right child; 2: virtual edge
			FLOAT sum, sum_resistance, tag, rPath, fPath;
		};
		nodetype *ds;
		
		void reset();
		FLOAT query(int x);
		void serere(int x, FLOAT value);
	};
	DS ds;
	
	struct otetype
	{
		int x,y;
		FLOAT resist, sum_r, flow;
	};
	otetype *ote;
	
	int *parent, *dfsorder;
	FLOAT *parente;
	
	FLOAT *str_sum;
	
	FLOAT *osum_rf, *tfPath, *tfdelta;
	
	Vec x,b;
	FLOAT err, bnorm;
	
	KOSZ();
	KOSZ(const GraphSP &g);
	void freeMemory();
	int sampleOTE();
	void toggle(int id);
	void setb(Vec _b);
	void convertFlowToSolution();
	void batchToggle(int order_num);
	SolverReturnValue solve(FLOAT tol, int maxbatch, int batchsize);
	SolverReturnValue solve(FLOAT tol, int maxbatch);
	SolverReturnValue solve(FLOAT tol);
	SolverReturnValue solve();
	operator AbstractSolver();
};

namespace KOSZHelper
{
	namespace HLD
	{
		struct btype
		{
			int which, lc, rc;
			FLOAT rtag;
		};
		
		vector< pair<int,FLOAT> > *e;

		struct data_type
		{
			FLOAT rPath, fPath;
			int which, id;
		};
		
		FLOAT *parente=new FLOAT[maxn];
		int *parent=new int[maxn], *size=new int[maxn], *weight=new int[maxn], *isHeavy=new int[maxn];
		int *isLeaf=new int[maxn], *s=new int[maxn], *tl=new int[maxn];
		btype *d2=new btype[maxn];
		
		int all;
		data_type *d=new data_type[maxn];
		
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
		
		pair<int,FLOAT> build(KOSZ::DS::nodetype *ds, int l, int r)
		{
			int root=find_mid(l,r); FLOAT s=0;
			int id=d[root].id; ds[id].tag=0; ds[id].sum=0;
			all++; tl[all]=id;
			if (root>l) 
			{
				pair<int,FLOAT> lc=build(ds, l,root-1); d2[id].lc=lc.first;
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
				pair<int,FLOAT> rc=build(ds, root+1,r); d2[id].rc=rc.first;
				ds[rc.first].parent=id; ds[rc.first].type=1; 
				ds[rc.first].rPath=d[root+1].rPath; ds[rc.first].fPath=0;
				d2[rc.first].which=d[root+1].which;
				s+=rc.second+d[root+1].rPath;
			}
			else	d2[id].rc=0;
			return make_pair(id,s+ds[id].sum_resistance);
		}

		void construct(KOSZ::DS::nodetype *ds, int m, int head)
		{
			s[0]=0;
			rep(i,1,m) s[i]=s[i-1]+weight[d[i].id];
			pair<int,FLOAT> rt=build(ds,1,m);
			ds[rt.first].parent=head; ds[rt.first].type=2; 
			ds[rt.first].rPath=d[1].rPath; ds[rt.first].fPath=0;
			d2[rt.first].which=d[1].which;
		}
		
		void hld(int n, vector< pair<int,FLOAT> > *_e, KOSZ::DS::nodetype *ds)
		{
			e=_e;
			memset(isHeavy,0,sizeof(int)*(n+1));
			memset(isLeaf,0,sizeof(int)*(n+1));
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
						d[k].rPath=parente[now]; 
						d[k].fPath=0; 
						d[k].which=now;
						d[k].id=now; 
						now=parent[now];
					}
					construct(ds,ci,z);
				}
		}
	}

	namespace DFS
	{
		vector< pair<int,FLOAT> > *e;
		int all;
		FLOAT *osum_r=new FLOAT[maxn];
		int *lis=new int[maxn];
		
		void dfs(int cur)
		{
			all++; lis[all]=cur;
			rept(it,e[cur])
			{
				osum_r[it->first]=osum_r[cur]+it->second;
				dfs(it->first);
			}
		}
		
		void setup(vector< pair<int,FLOAT> > *_e)
		{
			e=_e; osum_r[1]=0; all=0; dfs(1);
		}
	}
}

KOSZ::KOSZ()
{
	ds.ds=NULL;
	ote=NULL;
	parent=NULL;
	dfsorder=NULL;
	parente=NULL;
	str_sum=NULL;
	osum_rf=NULL;
	tfPath=NULL;
	tfdelta=NULL;
}

KOSZ::KOSZ(const GraphSP &g) 
{
	n=g.n; ds.n=n; ds.ds=new DS::nodetype[n+2]; parent=new int[n+2]; parente=new FLOAT[n+2];
	KOSZHelper::HLD::hld(n,g.e,ds.ds);
	rep(i,1,n) parent[i]=KOSZHelper::HLD::parent[i];
	rep(i,1,n) parente[i]=KOSZHelper::HLD::parente[i];
	
	osum_rf=new FLOAT[n+2]; dfsorder=new int[n+2]; 
	KOSZHelper::DFS::setup(g.e);
	rep(i,1,n) dfsorder[i]=KOSZHelper::DFS::lis[i];
	
	vector<FLOAT> ts=StretchCalculator::calculatePathResistance(g);
	m=g.o.size(); ote=new otetype[m+2]; str_sum=new FLOAT[m+2]; str_sum[0]=0;
	int all=0;
	rept(it,g.o) 
	{
		int x=get<0>(*it), y=get<1>(*it); FLOAT tres=ts[all]; all++; 
		ote[all].x=x; ote[all].y=y;
		ote[all].resist=get<2>(*it);
		ote[all].sum_r=ote[all].resist+tres;
		FLOAT stretch=tres/ote[all].resist;
		str_sum[all]=str_sum[all-1]+stretch+1;
	}
	
	tfPath=new FLOAT[n+2]; tfdelta=new FLOAT[n+2]; 
}

void KOSZ::freeMemory()
{
	delete[] ds.ds;
	delete[] ote;
	delete[] parent;
	delete[] dfsorder;
	delete[] parente;
	delete[] str_sum;
	delete[] osum_rf;
	delete[] tfPath;
	delete[] tfdelta;
}

int KOSZ::sampleOTE()
{
	FLOAT x=FLOAT(rand())/FLOAT(RAND_MAX)*str_sum[m];
	int t=lower_bound(str_sum+1,str_sum+m+1,x)-str_sum;
	return t-1;
}

void KOSZ::toggle(int id)
{
	int x=ote[id].x, y=ote[id].y;
	FLOAT sum_rf=ote[id].flow*ote[id].resist;
	FLOAT sum_r=ote[id].sum_r;
	sum_rf+=osum_rf[x]-osum_rf[y];
	sum_rf+=ds.query(x)-ds.query(y);
	FLOAT delta=sum_rf/sum_r;
	ote[id].flow-=delta;
	ds.serere(x,-delta); tfdelta[x]-=delta;
	ds.serere(y,delta); tfdelta[y]+=delta;
	//energy -= delta * sum_rf;
}

void KOSZ::setb(Vec _b)
{
	b=_b; assert(b.n==n);
	FLOAT s = 0;
	rep(i,0,n-1) s += b[i];
	s = s / FLOAT(n);
	rep(i,0,n-1) b[i] -= s;
	bnorm = b.norm();

	static FLOAT *sf=new FLOAT[maxn];
	rep(i,0,n) sf[i]=0;
	repd(kk,n,1)
	{
		int i=dfsorder[kk];
		sf[parent[i]]+=sf[i]-b[i-1];
		tfPath[i] = -(sf[i]-b[i-1]);
	}
	rep(i,1,m) ote[i].flow=0;
}
	
void KOSZ::convertFlowToSolution()
{
	static FLOAT *b1=new FLOAT[maxn];
	rep(i,1,n) b1[i]=0;
	x[0] = 0;
	rep(kk,1,n)
	{
		int i=dfsorder[kk];
		if (i==1) continue;
		FLOAT vv=tfPath[i]*parente[i];
		x[i-1] = x[parent[i]-1] + vv;
		b1[parent[i]]-=tfPath[i];
		b1[i]+=tfPath[i];
	}
	rep(i,1,m)
	{
		FLOAT fx=(x[ote[i].x-1] - x[ote[i].y-1]) / ote[i].resist;
		b1[ote[i].x]+=fx;
		b1[ote[i].y]-=fx;
	}
	
	FLOAT b1norm = 0;
	rep(i,0,n-1) b1norm+=(b1[i+1]-b[i])*(b1[i+1]-b[i]);
	b1norm=mysqrt(b1norm);
	err = b1norm/bnorm;
}
	
void KOSZ::batchToggle(int order_num)
{
	ds.reset();
	osum_rf[1]=0;
	rep(kk,1,n)
	{
		int i=dfsorder[kk];
		if (i>1) osum_rf[i]=osum_rf[parent[i]]+tfPath[i]*parente[i];
	}
	rep(i,0,n) tfdelta[i]=0;
	
	rep(i,1,order_num) toggle(sampleOTE()+1);
	
	repd(kk,n,1)
	{
		int i=dfsorder[kk];
		if (i>1) tfdelta[parent[i]]+=tfdelta[i];
	}
	rep(i,1,n) tfPath[i]+=tfdelta[i];
	
	convertFlowToSolution();
}
	
SolverReturnValue KOSZ::solve(FLOAT tol, int maxbatch, int batchsize)
{
	x=Vec(n);
	if (batchsize==-1) batchsize=4*n;
	SolverReturnValue ret;
	vector<FLOAT> &resvec=ret.resvec;
	convertFlowToSolution(); resvec.push_back(err);
	Vec bestx; bestx=x; FLOAT besterr=err; int whichit=0;
	int tt=0;
	while (1)
	{
		if (maxbatch!=-1 && tt>=maxbatch) break;
		batchToggle(batchsize);
		resvec.push_back(err); tt++;
		printf("KOSZ Iteration #%d relres = %.16lf\n",tt,err);
		if (maxbatch!=-1 && err<besterr) { besterr=err; bestx=x; whichit=tt; }
		if (err<tol) { ret.x=x; ret.flag=0; ret.relres=err; ret.iter=tt; return ret; };
	}
	ret.x=bestx; ret.flag=1; ret.relres=besterr; ret.iter=whichit; return ret;
}

SolverReturnValue KOSZ::solve(FLOAT tol, int maxbatch)
{
	return solve(tol,maxbatch,-1);
}

SolverReturnValue KOSZ::solve(FLOAT tol)
{
	return solve(tol,-1,-1);
}
	
SolverReturnValue KOSZ::solve()
{
	return solve(1e-6,-1,-1);
}

void KOSZ::DS::reset()
{
	rep(i,1,n)
	{
		ds[i].sum=0; ds[i].tag=0; ds[i].fPath=0;
	}
}

FLOAT KOSZ::DS::query(int x)
{
	FLOAT ans=0;
	while (x)
	{
		int flag=1; FLOAT rp=0;
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
		
void KOSZ::DS::serere(int x, FLOAT value)
{
	while (x)
	{
		int flag=1; FLOAT rp=0;
		while (1)
		{
			if (flag) 
			{
				ds[x].tag+=value;
				FLOAT tmp=ds[x].sum_resistance*value;
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

KOSZ::operator AbstractSolver()
{
	KOSZ t=*this;
	return AbstractSolver([t](const Vec &b, FLOAT tol, int maxit, const Vec &x0) mutable {
		t.setb(b); 
		return t.solve(tol,maxit);
	});
}

#endif
