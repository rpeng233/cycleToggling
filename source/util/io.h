/********************************************************************************
 * 
 * File IO Interface
 * 
 * Public API:
 * 	Mat IO::constructMatrixFromGraph(g):
 * 		GraphSP g: the graph
 * 		returns a Mat, the Laplacian matrix of graph g
 * 
 * 	GraphSP IO::convertLaplacianMatrixToGraphSP(const Mat &A)
 * 		Mat A: the Laplacian matrix
 * 		returns the graph corresponding to the Laplacian matrix
 * 		FIXME: 
 * 			Need to check input is indeed Laplacian matrix. 
 * 			Currently I'm not doing this bcuz of floating error issues.
 * 
 * 	GraphSP IO::readGraph(filename):
 * 		string filename: the input file name
 *		Reads input, calls the low stretch tree finder, and returns a GraphSP object.
 *   
 * 	GraphSP IO::readGraphSP(filename):
 * 		string filename: the input file name
 * 		Reads input, returns the GraphSP object.
 * 		It assumes that the first n-1 edges form a directed spanning tree rooted at 1.
 * 
 * 	Mat IO::readMML(filename):
 * 		string filename: the input file name
 * 		Reads the Matrix Market format input, and returns the matrix.
 *		It assumes that the input contains only the upper-triangular part of the matrix.
 * 
 * 	Mat IO::readMMA(filename):
 * 		string filename: the input file name
 * 		Reads the Matrix Market format input, and returns the matrix.
 *		It assumes that the input is the adjancy matrix of a graph, and contains only upper-triangular part.
 * 		It then constructs the Laplacian matrix corresponding to the graph.
 * 
 * 	Mat IO::readMMVec(filename):
 * 		string filename: the input file name
 * 		Reads the Matrix Market format input, and returns the vector.
 *		It assumes that the input is a dense matrix market format file, with a single column.
 *
 * 
 ********************************************************************************/

#ifndef __IO_H__
#define __IO_H__

#include "common.h"
#include "matrix.h"
#include "graph.h"
#include "treefinder.h"

namespace IO
{
	Mat constructMatrixFromGraph(const GraphSP &g)
	{
		int n=g.n;
		static FLOAT *s=new FLOAT[maxn];
#ifdef USE_MPFR
		rep(i,0,n) s[i]=0;
#else
		memset(s,0,sizeof(FLOAT)*(n+1));
#endif
		Mat A(n,n);
		rep(i,1,n)
			rept(it,g.e[i])
			{
				int x=i, y=it->first; FLOAT z=it->second;
				x--; y--; z=1.0/z;
				A.entryAddValue(x,y,-z); A.entryAddValue(y,x,-z);
				s[x]+=z; s[y]+=z;
			}
		rept(it,g.o)
		{
			int x=get<0>(*it), y=get<1>(*it); FLOAT z=get<2>(*it);
			x--; y--; z=1.0/z;
			A.entryAddValue(x,y,-z); A.entryAddValue(y,x,-z);
			s[x]+=z; s[y]+=z;
		}
		rep(i,0,n-1) A.entryAddValue(i,i,s[i]);
		A.sortup();
		return A;
	}
	
	GraphSP convertLaplacianMatrixToGraphSP(const Mat &A)	//won't report error if matrix is not symmetric
	{
		assert(A.n==A.m);
		Graph g(A.n);
		static FLOAT *s=new FLOAT[maxn];
		rept(it,A.values)
		{
			int x=it->x, y=it->y; FLOAT z=it->z;
			if (x!=y)
			{
				if (z>0) printf("Warning: Given matrix is not Laplacian, it has positive off-diagonals!\n");
				s[x]+=z; 
				z=-1.0/z;
				g.e[x+1].push_back(make_pair(y+1,z));
			}
			else
			{
				if (z<0) printf("Warning: Given matrix is not Laplacian, it has negative diagonals!\n");
				s[x]+=z;
			}
		}
		//rep(i,0,A.n-1) if (fabs(s[i])>1e-6) printf("Warning: Given matrix is not Laplacian, its row/col sum is not 0!\n");
		GraphSP g2=TreeFinder::findLowStretchTree(g);
		g.freeMemory();
		return g2;
	}
	
	GraphSP readGraphSP(string filename)
	{
		GraphSP g;
		FILE* ff = fopen(filename.c_str(), "r");
		if (!ff) { printf("File error.\n"); assert(0); }
		int n,m;
		fscanf(ff, "%d%d", &n, &m);
		char buff[100];
		fgets(buff,100,ff);	//skip the "order_num" parameter
		vector< pair<int,FLOAT> > *e=new vector< pair<int,FLOAT> >[n+1];
		vector< tuple<int,int,FLOAT> > o;
		rep(i,1,n) e[i].clear();
		rep(i,1,n-1)
		{
			int x,y; double z; fscanf(ff,"%d%d%lf",&x,&y,&z);
			e[x+1].push_back(make_pair(y+1,z));
		}
		o.clear();
		rep(i,1,m-(n-1))
		{
			int x,y; double z; fscanf(ff,"%d%d%lf",&x,&y,&z);
			o.push_back(make_tuple(x+1,y+1,z));
		}
		g.n=n; g.e=e; g.o.swap(o);
		fclose(ff);
		return g;
	}
	
	GraphSP readGraph(string filename)
	{
		FILE* ff = fopen(filename.c_str(), "r");
		if (!ff) { printf("File error.\n"); assert(0); }
		int n,m;
		fscanf(ff, "%d%d", &n, &m);
		char buff[100];
		fgets(buff,100,ff);	//skip the "order_num" parameter
		Graph h(n); 
		rep(i,1,n) h.e[i].clear();
		rep(i,1,m)
		{
			int x,y; double z; fscanf(ff,"%d%d%lf",&x,&y,&z);
			h.e[x+1].push_back(make_pair(y+1,z));
			h.e[y+1].push_back(make_pair(x+1,z));
		}
		GraphSP g=TreeFinder::findLowStretchTree(h);
		h.freeMemory();
		fclose(ff);
		return g;
	}

	//TOFIX Only works on symmetric mm right now
	Mat readMML(string filename)
	{
		FILE* ff = fopen(filename.c_str(), "r");
		if (!ff) { printf("File error.\n"); assert(0); }
		int c;
		char buff[100];
		
		do
		{
			fgets (buff, 100, ff);
			c = getc(ff);
			ungetc(c,ff);
		} while(c=='%');
		
		int n,temp,m;
		fscanf(ff, "%d%d%d", &n, &temp, &m);
				
		int x,y;
		double z;
		
		Mat A=Mat(n,n);
		rep(i,1,m)
		{
			fscanf(ff,"%d%d%lf",&x,&y,&z);
			x--; y--;
			A.entryAddValue(x,y,z);
			if (x!=y) A.entryAddValue(y,x,z);
		}
		
		fclose(ff);
			
		A.sortup();
		return A;
	}
	
	//TOFIX Only works on symmetric mm right now
	Mat readMMA(string filename)
	{
		FILE* ff = fopen(filename.c_str(), "r");
		if (!ff) { printf("File error.\n"); assert(0); }
		int c;
		char buff[100];
		
		do
		{
			fgets (buff, 100, ff);
			c = getc(ff);
			ungetc(c,ff);
		} while(c=='%');
		
		int n,temp,m;
		fscanf(ff, "%d%d%d", &n, &temp, &m);
		static FLOAT *s=new FLOAT[maxn];
#ifdef USE_MPFR
		rep(i,0,n) s[i]=0;
#else
		memset(s,0,sizeof(FLOAT)*(n+1));
#endif
		int x,y;
		double z;
		
		Mat A=Mat(n,n);
		rep(i,1,m)
		{
			fscanf(ff,"%d%d%lf",&x,&y,&z);
			x--; y--;

			A.entryAddValue(x,y,-z);
			A.entryAddValue(y,x,-z);
			s[x]+=z; s[y]+=z;
		}
		
		fclose(ff);
		
		rep(i,0,n-1) A.entryAddValue(i,i,s[i]);
		A.sortup();
		return A;
	}

        Vec readMMVec(string filename)
        {
		FILE* ff = fopen(filename.c_str(), "r");
		if (!ff) { printf("File error.\n"); assert(0); }
		int c;
		char buff[100];
		
		do
		{
			fgets (buff, 100, ff);
			c = getc(ff);
			ungetc(c,ff);
		} while(c=='%');
		
		int n,m;
		fscanf(ff, "%d%d", &n, &m);
		assert(m==1);

		Vec v(n);

		double z;
	     
		rep(i,0,n-1)
		{
			fscanf(ff,"%lf",&z);
			FLOAT zz=z;
			v[i]=zz;
		}
		
		fclose(ff);
		
		
		return v;
        }

	void saveMMVec(const Vec &A, string filename)
	{
		FILE* ff = fopen(filename.c_str(), "w");
		if (!ff) { printf("Write file error.\n"); assert(0); }
		fprintf(ff,"%%%%MatrixMarket matrix array real general\n");
		fprintf(ff,"%% Generated 01-Feb-2016\n");
		fprintf(ff,"%d %d\n",A.n,1);
		rep(i,0,A.n-1) fprintf(ff,"%.16lf\n",printFloat(A[i]));
		fclose(ff);
	}
	
	
	void saveMM(const Mat &A, string filename)
	{
		A.sortup();
		FILE* ff = fopen(filename.c_str(), "w");
		if (!ff) { printf("Write file error.\n"); assert(0); }
		fprintf(ff,"%%%%MatrixMarket matrix coordinate real symmetric\n");
		fprintf(ff,"%% Generated 01-Feb-2016\n");
		int cnt=0;
		rept(it,A.values) if (it->x>=it->y) cnt++;
		fprintf(ff,"%d %d %d\n",A.n,A.m,cnt);
		rept(it,A.values) if (it->x>=it->y) fprintf(ff,"%d %d %.16lf\n",it->x+1,it->y+1,printFloat(it->z));
		fclose(ff);
	}
	
	void saveMMnonsym(const Mat &A, string filename)
	{
		A.sortup();
		FILE* ff = fopen(filename.c_str(), "w");
		if (!ff) { printf("Write file error.\n"); assert(0); }
		fprintf(ff,"%%%%MatrixMarket matrix coordinate real general\n");
		fprintf(ff,"%% Generated 01-Feb-2016\n");
		int cnt=0;
		rept(it,A.values) cnt++;
		fprintf(ff,"%d %d %d\n",A.n,A.m,cnt);
		rept(it,A.values) fprintf(ff,"%d %d %.16lf\n",it->x+1,it->y+1,printFloat(it->z));
		fclose(ff);
	}
		
	//assumes non-symmetric!
	Mat readMMnonsym(string filename)
	{
		FILE* ff = fopen(filename.c_str(), "r");
		if (!ff) { printf("File error.\n"); assert(0); }
		int c;
		char buff[100];
		
		do
		{
			fgets (buff, 100, ff);
			c = getc(ff);
			ungetc(c,ff);
		} while(c=='%');
		
		int n,m,nnz;
		fscanf(ff, "%d%d%d", &n, &m, &nnz);
				
		int x,y;
		double z;
		
		Mat A=Mat(n,m);
		rep(i,1,nnz)
		{
			fscanf(ff,"%d%d%lf",&x,&y,&z);
			x--; y--;
			A.entryAddValue(x,y,z);
		}
		
		fclose(ff);
			
		A.sortup();
		return A;
	}
}
	  
	

#endif
