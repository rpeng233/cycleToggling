/********************************************************************************
 * 
 * Low Stretch Spanning Tree utility
 * 
 * Public API:
 * 	vector<FLOAT> StretchCalculator::calculatePathResistance(const GraphSP &g):
 * 		Input: g, the graph.
 * 		Returns: 
 * 			The tree-path resistance corresponding to off-tree edges of g, in the same order as g.o.
 *			Stretch of an off-tree edge is by definition (tree-path resistance)/(off-tree edge resistance).
 *	 
 * 	FLOAT StretchCalculator::calculateTotalStretch(const GraphSP &g):
 * 		Input: g, the graph.
 * 		Returns: the total stretch of g.
 * 
 * 	GraphSP TreeFinder::findLowStretchTree(const Graph &g):
 * 		Input: g, the graph.
 * 		Returns: 
 * 			same graph, but with a spanning tree.
 * 			currently it tries Dijkstra tree and MST, and picks the better one.
 * 
 * 	GraphSP TreeFinder::findLowStretchTree(const GraphSP &g):
 * 		Input: g, the graph.
 * 		Returns: 
 * 			same graph, but with a spanning tree.
 * 			It ignores the current spanning tree in g, and finds a new one.
 * 
 ********************************************************************************/

#ifndef __TREEFINDER_H__
#define __TREEFINDER_H__

#include "common.h"
#include "graph.h"

namespace TarjanLCA
{
	vector< pair<int,FLOAT> > *e;
	vector< pair<int,int> > *q=new vector< pair<int,int> >[maxn];
	vector<int> ret;
	int *p=new int[maxn], *vis=new int[maxn];
	FLOAT *osum_r=new FLOAT[maxn];
	
	int find(int x)
	{
		if (p[x]==x) return x;
		p[x]=find(p[x]);
		return p[x];
	}
		
	void dfs(int cur)
	{
		rept(it,e[cur])
		{
			osum_r[it->first]=osum_r[cur]+it->second;
			dfs(it->first);
			p[find(it->first)]=find(cur);
		}
		vis[cur]=1;
		rept(it,q[cur])
			if (vis[it->first])
				ret[it->second]=find(it->first);
	}
			
	vector<int> solve(int n, vector< pair<int, FLOAT> > *_e, const vector< tuple<int,int,FLOAT> > &o)
	{
		e=_e;
		int all=o.size();
		rep(i,1,n) q[i].clear();
		rep(i,0,all-1)
		{
			int x=get<0>(o[i]), y=get<1>(o[i]);
			q[x].push_back(make_pair(y,i));
			q[y].push_back(make_pair(x,i));
		}
		rep(i,1,n) vis[i]=0;
		rep(i,1,n) p[i]=i;
		ret.resize(o.size());
		osum_r[1]=0; dfs(1);
		vector<int> oret=ret;
		return oret;
	}
}

namespace StretchCalculator
{
	vector<FLOAT> calculatePathResistance(const GraphSP &g)
	{
		vector<int> lca=TarjanLCA::solve(g.n,g.e,g.o);
		vector<FLOAT> ret; int all=0;
		rept(it,g.o) 
		{
			int x=get<0>(*it), y=get<1>(*it), z=lca[all]; all++;
			ret.push_back(TarjanLCA::osum_r[x]+TarjanLCA::osum_r[y]-TarjanLCA::osum_r[z]*2);
		}
		return ret;
	}
	
	FLOAT calculateTotalStretch(const GraphSP &g)
	{
		vector<FLOAT> t=calculatePathResistance(g);
		FLOAT ret=0; int all=0;
		rept(it,g.o) { ret+=t[all]/get<2>(*it); all++; }
		return ret;
	}
}

namespace DijkstraTreeFinder
{
	int *place=new int[maxn], *prev=new int[maxn];
	pair<int,FLOAT> *a=new pair<int,FLOAT>[maxn];
	FLOAT *dis=new FLOAT[maxn], *prevd=new FLOAT[maxn];
	int hps;
	
	void heapdown(int i)
	{
		while (1)
		{
			int mini;
			if (i*2<=hps && a[i*2].second<a[i].second) mini=i*2; else mini=i;
			if (i*2+1<=hps && a[i*2+1].second<a[mini].second) mini=i*2+1;
			if (i!=mini)
			{
				swap(a[i],a[mini]);
				place[a[i].first]=i; place[a[mini].first]=mini;
				i=mini;
			}
			else	break;
		}
	}
	
	void heapup(int i)
	{
		while (i>1 && a[i/2].second>a[i].second)
		{
			swap(a[i],a[i/2]);
			place[a[i].first]=i; place[a[i/2].first]=i/2;
			i/=2;
		}
	}
	
	void insert(pair<int,FLOAT> key)
	{
		hps++; a[hps]=key; place[a[hps].first]=hps;
		heapup(hps);
	}
	
	int extractMin()
	{
		int ret=a[1].first; 
		a[1]=a[hps]; place[a[1].first]=1; hps--;
		heapdown(1);
		return ret;
	}
	
	void dfs(GraphSP &h, int cur)
	{
		rept(it,h.e[cur])
		{
			rept(it2,h.e[it->first])
				if (it2->first==cur)
				{
					h.e[it->first].erase(it2); break;
				}
			dfs(h,it->first);
		}
	}
	
	GraphSP solve(const Graph &g)
	{
		int source=rand()%g.n+1;
		hps=1; a[1]=make_pair(source,0); 
		rep(i,1,g.n) place[i]=-1; place[source]=1;
		rep(i,1,g.n) dis[i]=1e100; dis[source]=0;
		prev[source]=-1;
		while (hps)
		{
			int x=extractMin(); place[x]=-2;
			rept(it,g.e[x])
				if (place[it->first]!=-2)
				{
					FLOAT w=it->second;
					if (dis[x]+w<dis[it->first])
					{
						dis[it->first]=dis[x]+w;
						prev[it->first]=x; prevd[it->first]=it->second;
						if (place[it->first]!=-1)
						{
							a[place[it->first]].second=dis[it->first];
							heapup(place[it->first]);
						}
						else	insert(make_pair(it->first,dis[it->first]));
					}
				}
		}
		rep(i,1,g.n) assert(dis[i]<1e99);
		GraphSP h(g.n);
		rep(i,1,g.n)
			if (prev[i]!=-1)
			{
				h.e[i].push_back(make_pair(prev[i],prevd[i]));
				h.e[prev[i]].push_back(make_pair(i,prevd[i]));
			}
		dfs(h,1);
		rep(i,1,g.n)
			rept(it,g.e[i])
				if (prev[i]!=it->first && prev[it->first]!=i)
					if (i<it->first)
						h.o.push_back(make_tuple(i,it->first,it->second));
		return h;
	}
}

namespace MSTTreeFinder
{
	int *p=new int[maxn], *prev=new int[maxn];
	
	int find(int x)
	{
		if (p[x]==x) return x;
		p[x]=find(p[x]);
		return p[x];
	}

	struct atype
	{
		int u,v;
		FLOAT w;
		int r;
	};

	int cmp(const atype &a, const atype &b)
	{
		if (a.w!=b.w) return a.w<b.w;
		if (a.r!=b.r) return a.r>b.r;
		if (a.u!=b.u) return a.u<b.u;
		return a.v<b.v;
	}
	
	void dfs(GraphSP &h, int cur)
	{
		rept(it,h.e[cur])
		{
			rept(it2,h.e[it->first])
				if (it2->first==cur)
				{
					h.e[it->first].erase(it2); break;
				}
			prev[it->first]=cur;
			dfs(h,it->first);
		}
	}
	
	GraphSP solve(const Graph &g)
	{
		int n=g.n;
		vector<atype> lis;
		rep(i,1,n)
			rept(it,g.e[i])
			{
				atype t; t.u=i; t.v=it->first; t.w=it->second; t.r=rand();
				lis.push_back(t);
			}
		sort(lis.begin(),lis.end(),cmp);
		
		GraphSP h(n);
		rep(i,1,n) p[i]=i;
		int cnt=0;
		rept(it,lis)
			if (find(it->u)!=find(it->v))
			{
				//printf("%d %d\n",it->u,it->v);
				p[find(it->u)]=find(it->v); 
				h.e[it->u].push_back(make_pair(it->v,it->w));
				h.e[it->v].push_back(make_pair(it->u,it->w));
				cnt++; if (cnt==n-1) break;
			}	
		assert(cnt==n-1);
		
		prev[1]=-1; dfs(h,1);
		rep(i,1,g.n)
			rept(it,g.e[i])
				if (prev[i]!=it->first && prev[it->first]!=i)
					if (i<it->first)
						h.o.push_back(make_tuple(i,it->first,it->second));
		
		return h;
	}
}

namespace TreeFinder
{
	GraphSP findLowStretchTree(const Graph &g)
	{
		static FLOAT *vis=new FLOAT[maxn];
		rep(i,1,g.n) vis[i]=0;
		Graph g2(g.n);
		rep(i,1,g.n)
		{
			rept(it,g.e[i]) vis[it->first]+=1.0/it->second;
			rept(it,g.e[i])
				if (vis[it->first]!=0)
				{
					g2.e[i].push_back(make_pair(it->first,1.0/vis[it->first]));
					vis[it->first]=0;
				}
		}
		GraphSP h1=DijkstraTreeFinder::solve(g2);
		GraphSP h2=MSTTreeFinder::solve(g2);
		g2.freeMemory();
		FLOAT st1=StretchCalculator::calculateTotalStretch(h1);
		FLOAT st2=StretchCalculator::calculateTotalStretch(h2);
		if (st1<st2)
		{
			h2.freeMemory();
			//printf("findLowStretchTree (Dijkstra): average stretch = %.16lf\n",(double)(st1/h1.o.size()));
			return h1;
		}
		else
		{
			h1.freeMemory();
			//printf("findLowStretchTree (MST): average stretch = %.16lf\n",(double)(st2/h2.o.size()));
			return h2;
		}
	}
	
	GraphSP findLowStretchTree(const GraphSP &g)
	{
		Graph g2=Graph(g.n);
		rep(i,1,g.n) 
			rept(it,g.e[i])
			{
				g2.e[i].push_back(*it);
				g2.e[it->first].push_back(make_pair(i,it->second));
			}
		rept(it,g.o)
		{
			int x=get<0>(*it), y=get<1>(*it); FLOAT z=get<2>(*it);
			g2.e[x].push_back(make_pair(y,z));
			g2.e[y].push_back(make_pair(x,z));
		}
		GraphSP h=findLowStretchTree(g2);
		g2.freeMemory();
		return h;
	}
}

#endif
