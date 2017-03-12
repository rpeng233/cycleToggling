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

GraphSP UltraSparsifier(const GraphSP &g, double kappa, double c, double lgn)
{
	int n=g.n;
	GraphSP h(n); 
	rep(i,1,n) 
		rept(it,g.e[i])
			h.e[i].push_back(make_pair(it->first,it->second/kappa));
	
	//FLOAT sm=StretchCalculator::calculateTotalStretch(g)*lgn*c/kappa;
	int sm=int(n/5);
	
	vector<FLOAT> str=StretchCalculator::calculatePathResistance(g);
	FLOAT *s=new FLOAT[str.size()+2];
	s[0]=0; int all=0;
	rept(it,g.o)
	{
		FLOAT est=str[all]/get<2>(*it); //if (est>kappa) est=kappa;
		all++; s[all]=s[all-1]+est;
	}
	
	set<int> ss;
	int m=str.size();
	FLOAT i=0;
	while (i<sm)
	{
		i=i+1.0;
		FLOAT x=FLOAT(rand())/FLOAT(RAND_MAX)*s[m];
		int t=lower_bound(s+1,s+m+1,x)-s;
		t--;
		FLOAT nr=c*lgn*str[t]/kappa;
		//if (ss.find(t)!=ss.end()) continue;
		ss.insert(t);
		h.o.push_back(make_tuple(get<0>(g.o[t]),get<1>(g.o[t]),nr));
	}
	delete[] s;
	return h;
	//return TreeFinder::findLowStretchTree(h);
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

	