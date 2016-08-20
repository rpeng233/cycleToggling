/********************************************************************************
 * 
 * Ultra Sparsifier
 * 
 * Public API:
 * 	UltraSparsifier(const GraphSP &g, int m):
 * 		Inputs:
 * 			g: the input graph
 * 			m: the off-tree edges to be sampled
 * 		
 * 		Returns:
 * 			the result graph
 * 
 * Note:
 * 	Currently the sparsifier just samples m edge according to stretch, 
 * 	and remove any duplicated edges. So the result might contain less than m edges.
 * 
 *******************************************************************************/

#ifndef __ULTRASPARSIFIER_H__
#define __ULTRASPARSIFIER_H__

#include "common.h"
#include "matrix.h"
#include "io.h"
#include "kosz.h"
#include "treefinder.h"

GraphSP UltraSparsifier(const GraphSP &g, double scale, int sm)
{
	int n=g.n;
	GraphSP h(n); 
	rep(i,1,n) 
		rept(it,g.e[i])
			h.e[i].push_back(make_pair(it->first,it->second*scale));
	
	vector<FLOAT> str=StretchCalculator::calculatePathResistance(g);
	int all=0; FLOAT sumStr=0;
	rept(it,g.o)
	{
		str[all]=str[all]/get<2>(*it); sumStr+=str[all];
		all++; 
	}
	
	all=0;
	rept(it,g.o)
	{
		FLOAT p_e=sm*str[all]/sumStr;
		FLOAT z=FLOAT(rand())/FLOAT(RAND_MAX);
		if (z<p_e)  h.o.push_back(make_tuple(get<0>(*it),get<1>(*it),get<2>(*it)*p_e));
		all++;
	}
	return h;
	return TreeFinder::findLowStretchTree(h);
}

GraphSP UltraSparsifier(const GraphSP &g, int sm)
{
	KOSZ z=KOSZ(g);	//just to sample off-tree edges....
	vector< tuple<int,int,FLOAT> > no;
	set<int> vis;
	rep(i,1,sm)
	{
		int t=z.sampleOTE();
		if (vis.find(t)!=vis.end()) continue;
		vis.insert(t);
		no.push_back(g.o[t]);
	}
	z.freeMemory(); 
	int n=g.n;
	GraphSP h(n); rep(i,1,n) h.e[i]=g.e[i]; h.o.swap(no);
	return h;
}

#endif

	