/********************************************************************************
 * 
 * Graph Contraction Solver
 * 
 * On a connected graph G with n vertices and m edges, 
 * G can be contracted to a graph H with no more than 4(m-n+1) vertices and 5(m-n+1) edges
 * A Laplacian solver for H can be converted to a solver for G in linear time.
 * Graph Contraction Solver takes a graph G and a solver for H, and constructs the solver for G.
 * 
 * Public API:
 * 	ContractSolver(const GraphSP &g, GraphSP &h):
 *		Inputs:
 *			GraphSP g, the input graph
 * 		
 * 		Outputs:
 * 			The contracted graph of g is stored in h
 * 			Returns a function f, which takes an AbstractSolver as parameter, and returns an AbstractSolver
 * 			If f is fed with an AbstractSolver for h, it returns an AbstractSolver for g.
 * 
 * 		Example:
 * 			auto f=ContractSolver(g,h);	
 * 			AbstractSolver S=KOSZ(h);	//construct KOSZ solver for graph h
 * 			AbstractSolver T=f(S);		//now T is a solver for graph g
 * 
 * 	ContractSolver(const GraphSP &g, const function<AbstractSolver(const GraphSP&)> &f):
 * 		Inputs:
 *			GraphSP g, the input graph
 * 			f: a function which, when given a graph as input, returns an AbstractSolver for that graph.
 * 
 * 		Outputs:
 * 			the AbstractSolver for g.
 * 
 * 		Example:
 * 			AbstractSolver S=ContractSolver(g,[](const GraphSP &h) {
 * 				return (AbstractSolver)KOSZ(h);
 *			});
 * 
 ********************************************************************************/

#ifndef __CONTRACTSOLVER_H__
#define __CONTRACTSOLVER_H__

#include "common.h"
#include "matrix.h"
#include "abstractsolver.h"
#include "graph.h"
#include "treefinder.h"

namespace ContractSolverHelper
{
	void dfsTree(const GraphSP &g, int cur, int *mark, vector<int> &res)
	{
		rept(it,g.e[cur])
			if (mark[it->first]==0)
			{
				res.push_back(it->first);
				dfsTree(g,it->first,mark,res);
			}
	}

	vector<int> findCriticalNodes(const GraphSP &g)
	{
		int n=g.n;
		static int *isCritical=new int[maxn];
		memset(isCritical,0,sizeof(int)*(n+1));
		rept(it,g.o)
		{
			int x=get<0>(*it), y=get<1>(*it);
			isCritical[x]=1; isCritical[y]=1;
		}
		
		static int *q=new int[maxn];
		int head=1, tail=2; q[head]=1;
		while (head<tail)
		{
			int x=q[head]; head++;
			rept(it,g.e[x]) 
			{ 
				q[tail]=it->first; tail++; 
			}
		}
		static int *f=new int[maxn];
		memset(f,0,sizeof(int)*(n+1));
		repd(kk,n,1)
		{
			int i=q[kk];
			//flag: 0 not exist, >0 critical child id. real node has flag==i
			int flag=0;
			if (isCritical[i]) flag=i;
			if (flag==0)
			{
				int cnt=0;
				rept(it,g.e[i])
					if (f[it->first]>0) { cnt++; flag=f[it->first]; }
				if (cnt>=2 || i==1) flag=i;
			}
			f[i]=flag;
		}
		
		vector<int> ret;
		rep(i,1,n) if (f[i]==i) ret.push_back(i);
		return ret;
	}

	void contract(const GraphSP &g, GraphSP &h, vector<int> &lis)
	{
		int n=g.n;
		lis=findCriticalNodes(g);
		int all=0;
		static int *which=new int[maxn];
		memset(which,0,sizeof(int)*(n+1));
		rept(it,lis) { all++; which[*it]=all; }
		static int *p=new int[maxn]; static FLOAT *pw=new FLOAT[maxn];
		rep(i,1,n) rept(it,g.e[i]) { p[it->first]=i; pw[it->first]=it->second; }
		h=GraphSP(all);
		rept(it,lis) 
		{
			int i=*it; FLOAT w=0;
			if (i==1) continue;
			while (which[p[i]]==0) { w+=pw[i]; i=p[i]; }
			w+=pw[i];
			h.e[which[p[i]]].push_back(make_pair(which[*it],w));
		}
		rept(it,g.o) h.o.push_back(make_tuple(which[get<0>(*it)], which[get<1>(*it)], get<2>(*it)));
	}
}

function<AbstractSolver(const AbstractSolver&)> ContractSolver(const GraphSP &g, GraphSP &hs)
{
	int n=g.n; GraphSP h; vector<int> tmp;
	ContractSolverHelper::contract(g,h,tmp);
	int *mp=new int[h.n+1]; mp[0]=0; rep(i,1,h.n) mp[i]=tmp[i-1];
	int *which=new int[n+1]; memset(which,0,sizeof(int)*(n+1));
	rep(i,1,h.n) which[mp[i]]=i;
	
	static int *ph=new int[maxn];
	int *p=new int[n+1];
	FLOAT *pw=new FLOAT[n+1];
	FLOAT *phw=new FLOAT[h.n+1];
	rep(i,1,n) rept(it,g.e[i]) p[it->first]=i;
	rep(i,1,n) rept(it,g.e[i]) pw[it->first]=it->second;
	rep(i,1,h.n) rept(it,h.e[i]) ph[it->first]=i;
	rep(i,1,h.n) rept(it,h.e[i]) phw[it->first]=it->second; phw[1]=0;
	
	vector< pair<int,FLOAT> > *vlis=new vector< pair<int,FLOAT> >[h.n+1];
	static int *mark=new int[maxn];
	rep(i,1,n) mark[i]=0;
	rep(kk,1,h.n)
	{
		int i=mp[kk]; mark[i]=1; FLOAT x=0;
		if (i==1) continue;
		while (i!=mp[ph[kk]])
		{
			x+=pw[i]; i=p[i]; mark[i]=1;
			vlis[kk].push_back(make_pair(i,x));
		}
		reverse(vlis[kk].begin(),vlis[kk].end());
		int ss=int(vlis[kk].size())-1;
		rep(i,1,ss) vlis[kk][i].second/=vlis[kk][0].second;
	}
	
	vector<int> *leafList=new vector<int>[n+1];
	rep(i,1,n)
	{
		leafList[i].clear();
		if (mark[i]) ContractSolverHelper::dfsTree(g,i,mark,leafList[i]);
	}
	
	int hn=h.n;
	FLOAT *ps=new FLOAT[h.n+1];
	int *dfsorder=new int[h.n+1];
	int head=1, tail=2; dfsorder[head]=1;
	while (head<tail)
	{
		int x=dfsorder[head]; head++;
		rept(it,h.e[x]) { dfsorder[tail]=it->first; tail++; }
	}
	
	hs=TreeFinder::findLowStretchTree(h);
	//Mat mH=IO::constructMatrixFromGraph(h);
	h.freeMemory();
	
	return [n,hn,which,mp,p,dfsorder,phw,pw,ps,vlis,leafList](const AbstractSolver &S) {
		return AbstractSolver([=](const Vec &_b, FLOAT tol, int maxit) {
			assert(_b.n==n); Vec b=_b; 
			Vec nb(hn);
			rep(i,1,hn) ps[i]=0;
			rep(i,1,n)
				repd(j,int(leafList[i].size())-1,0)
					b[p[leafList[i][j]]-1]+=b[leafList[i][j]-1];

			rep(i,1,hn) 
			{
				nb[i-1]+=b[mp[i]-1];
				if (i==1) continue;
				int ss=int(vlis[i].size())-1;
				FLOAT x1=0, x2=0;
				rep(j,1,ss) 
				{
					x1+=b[vlis[i][j].first-1]*vlis[i][j].second;
					x2+=b[vlis[i][j].first-1];
				}
				ps[i]=x1; nb[which[vlis[i][0].first]-1]+=x1; nb[i-1]+=x2-x1;
			}
			
			//Vec x0(hn);
			//rep(i,1,hn) x0[i-1]=_x0[mp[i]-1];
			
			SolverReturnValue result=S.solve(nb,tol,maxit);	
			
			Vec &x=result.x;
			Vec ret(n); 
			rep(j,1,hn)
			{
				int kk=dfsorder[j];
				if (mp[kk]==1) continue;
				FLOAT flow=(x[which[vlis[kk][0].first]-1]-x[kk-1])/phw[kk]-ps[kk];
				FLOAT st=ret[vlis[kk][0].first-1];
				rep(i,1,int(vlis[kk].size())-1)
				{
					st+=pw[vlis[kk][i].first]*flow;
					ret[vlis[kk][i].first-1]=st;
					flow+=b[vlis[kk][i].first-1];
				}
				ret[mp[kk]-1]=st+pw[mp[kk]]*flow;
			}
			rep(i,1,n)
				rept(it,leafList[i])
				{
					int k=*it;
					ret[k-1]=ret[p[k]-1]-b[k-1]*pw[k];
				}
			result.x=ret*FLOAT(-1);
			return result;
		});
	};
}

AbstractSolver ContractSolver(const GraphSP &g, const function<AbstractSolver(const GraphSP&)> &Sv)
{
	GraphSP h; auto fn=ContractSolver(g,h);
	return fn(Sv(h));
}

#endif