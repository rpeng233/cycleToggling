/********************************************************************************
 * 
 * Graph Definitions
 * 
 * API Interface:
 * 	struct Graph: definition of a graph (1-INDEXED).
 * 		int .n: the size of graph
 * 		vector< pair<int,FLOAT> > .e[i]: the neighbors of vertex i, 1<=i<=n
 * 		.freeMemory(): destroys graph and frees memory
 * 
 * 	struct GraphSP: definition of a connected graph with a spanning tree (1-INDEXED).
 * 		int .n: the size of graph
 * 		vector< pair<int,FLOAT> > .e[i]: the children of vertex i in spanning tree, 1<=i<=n
 * 		                                  the tree is DIRECTED, and 1 is always the root
 * 		vector< tuple<int,int,FLOAT> > .o: the off-tree edges
 * 		.freeMemory(): destroys graph and frees memory
 * 
 * NOTE:
 * 	The weight in Graph and GraphSP objects are always RESISTANCE, not WEIGHT!!!
 * 
 * 	Graph and GraphSP behave like objects in python or javascript. That is, if you do 
 * 		Graph A; Graph B=A; 
 * 	Then A and B will be actually pointing to the same object, 
 * 	i.e. modifying B will result in A being modified as well!
 * 	
 * 	However, C++ does not have an automatic garbage collection system, 
 * 	so you have to run .freeMemory() to free memory of graphs that are no longer needed.
 *
 * 	In Graph, every edge appears twice, in both direction.
 * 	In GraphSP, however, every edge appears in only one direction.
 * 
 ********************************************************************************/

#ifndef __GRAPH_H__
#define __GRAPH_H__

#include "common.h"

struct Graph
{
	int n;
	vector< pair<int, FLOAT> > *e;
	Graph(): n(0), e(NULL) {}
	Graph(int n): n(n), e(new vector< pair<int,FLOAT> >[n+1]) {}
	Graph(const Graph &b): n(b.n), e(b.e) {}
	Graph &operator =(const Graph&& b) { n=b.n; e=b.e; return (*this); }
	void freeMemory() const { delete[] e; }
};

struct GraphSP
{
	int n;
	vector< pair<int, FLOAT> > *e;
	vector< tuple<int,int,FLOAT> > &o;
	GraphSP(): n(0), e(NULL), o(*new vector< tuple<int,int,FLOAT> >) {}
	GraphSP(int n): n(n), e(new vector< pair<int,FLOAT> >[n+1]), o(*new vector< tuple<int,int,FLOAT> >) {}
	GraphSP(const GraphSP &b): n(b.n), e(b.e), o(b.o) {}
	GraphSP &operator =(const GraphSP&& b) { n=b.n; e=b.e; o.swap(b.o); return (*this); }
	void freeMemory() const { delete &o; delete[] e; }
};

#endif