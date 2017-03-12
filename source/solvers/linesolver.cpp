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

/**
 * please adjust bitn so that 2^bitn > (n+10)/block (do NOT ignore +10)
 * maxl > m+10
 */

#define block 20

int bitn, offset;

#define maxm 3000010
#define maxl 4000010

bool init_run = true;

char msg[10000];

struct data_type
{
	double fPath, rPath, ts;
};
data_type d[maxl];

namespace DS
{
	struct atype
	{
		double flow, sum_resistance, sum;
	};
	
	atype a[maxn];
	
	void init(data_type *z, int n)
	{
		memset(a,0,offset*2*sizeof(atype));
		rep(i,1,n/block) 
		{
			double s=0;
			rep(j,(i-1)*block+1,i*block) s+=z[j].rPath;
			a[offset+i].sum_resistance=s;
		}
		repd(i,offset-1,1) a[i].sum_resistance=a[i*2].sum_resistance+a[i*2+1].sum_resistance;
	}
	
	void apply(int i, double value)
	{
		a[i].flow+=value;
		a[i].sum+=value*a[i].sum_resistance;
	}
	
	void serere(int l, int r, double value, double _lp, double _rp)
	{
		l--; r++; l+=offset; r+=offset;
		double lp=_lp, rp=_rp;
		while (1)
		{
			a[l].sum+=lp*value; a[r].sum+=rp*value;
			if (r-l==1) break;
			if (l%2==0) apply(l+1,value), lp+=a[l+1].sum_resistance; 
			if (r%2==1) apply(r-1,value), rp+=a[r-1].sum_resistance; 
			l/=2; r/=2;
		}
		l/=2; while (l>0) a[l].sum+=lp*value, l/=2;
		r/=2; while (r>0) a[r].sum+=rp*value, r/=2;
	}
	
	void serere_single(int x, double s)
	{
		x+=offset;
		while (x>0)
		{
			a[x].sum+=s; x/=2;
		}
	}
	
	double query_single(int x)
	{
		double res=0;
		x+=offset;
		while (x>0)
		{
			res+=a[x].flow; x/=2;
		}
		return res;
	}
	
	double query(int l, int r)
	{
		l--; r++; l+=offset; r+=offset;
		double ans=0, lp=0, rp=0;
		while (1)
		{
			ans+=lp*a[l].flow+rp*a[r].flow;
			if (r-l==1) break;
			if (l%2==0) ans+=a[l+1].sum, lp+=a[l+1].sum_resistance;
			if (r%2==1) ans+=a[r-1].sum, rp+=a[r-1].sum_resistance; 
			l/=2; r/=2;
		}
		l/=2; while (l>0) ans+=lp*a[l].flow, l/=2;
		r/=2; while (r>0) ans+=rp*a[r].flow, r/=2;
		return ans;
	}
		
	void dump(double *res, int n)	//this destroys DS!
	{
		rep(i,1,offset-1)
		{
			a[i*2].flow+=a[i].flow;
			a[i*2+1].flow+=a[i].flow;
		}
		rep(i,1,n/block) 
			rep(j,(i-1)*block+1,i*block)
				res[j]=a[offset+i].flow;
	}
}

double resistance[maxl], f[maxl];
double osum_r[maxl], osum_rf[maxl], ans[maxl];
int order[maxl], e[maxl][2];

void serere(int l, int r, double value)
{
	double lp=0, rp=0;
	int sl=(l-2)/block*block+block+1;
	if (l==1) sl=1;
	int sr=r/block*block;
	if (sl<sr)
	{
		while (l<sl)
		{
			d[l].fPath+=value; d[l].ts+=value*d[l].rPath; lp+=d[l].rPath; l++;
		}
		while (r>sr)
		{
			d[r].fPath+=value; d[r].ts+=value*d[r].rPath; rp+=d[r].rPath; r--;
		}
		DS::serere((l-1)/block+1,r/block,value,lp,rp);
	}
	else  if (sl==sr+1)
	{
		rep(i,l,sr)
		{
			double t=value*d[i].rPath;
			d[i].fPath+=value; d[i].ts+=t; lp+=t;
		}
		rep(i,sl,r)
		{
			double t=value*d[i].rPath;
			d[i].fPath+=value; d[i].ts+=t; rp+=t;
		}
		DS::serere_single(sr/block,lp);
		DS::serere_single(sr/block+1,rp);
	}
	else
	{
		rep(i,l,r)
		{
			double t=value*d[i].rPath;
			d[i].fPath+=value; d[i].ts+=t; lp+=t;
		}
		DS::serere_single(sl/block,lp);
	}
}

double query(int l, int r)
{
	double ret=0, lp=0, rp=0;
	int sl=(l-2)/block*block+block+1;
	if (l==1) sl=1;
	int sr=r/block*block;
	if (sl<=sr)
	{
		while (l<sl)
		{
			ret+=d[l].ts; lp+=d[l].rPath; l++;
		}
		while (r>sr)
		{
			ret+=d[r].ts; rp+=d[r].rPath; r--;
		}
		ret+=DS::query((l-1)/block+1,r/block);
		if (lp!=0) ret+=lp*DS::query_single(l/block);
		if (rp!=0) ret+=rp*DS::query_single(r/block+1);
		return ret;
	}
	else  if (sl==sr+1)
	{
		rep(i,l,sr)
		{
			ret+=d[i].ts; lp+=d[i].rPath;
		}
		rep(i,sl,r)
		{
			ret+=d[i].ts; rp+=d[i].rPath;
		}
		ret+=lp*DS::query_single(sl/block);
		ret+=rp*DS::query_single(sl/block+1);
		return ret;
	}
	else
	{
		rep(i,l,r) { ret+=d[i].ts; lp+=d[i].rPath; }
		ret+=lp*DS::query_single(sl/block);
		return ret;
	}
}

double energy;

void toggle(int id)
{
	int l=e[id][0]+1, r=e[id][1];
	double sum_rf=-f[id]*resistance[id];
	double sum_r=resistance[id];
	sum_rf+=osum_rf[r]-osum_rf[l-1];
	sum_r+=osum_r[r]-osum_r[l-1];
	sum_rf+=query(l,r);
	//printf("%d %d\n%.3lf %.3lf\n",l,r,sum_rf,sum_r);
	double delta=sum_rf/sum_r;
	f[id]+=delta;
	serere(l,r,-delta);
	energy -= delta * sum_rf;
}

double x[maxl], b[maxl], b1[maxl];

int n,m, order_num, tt=0;
double bnorm;

double sqr(double x) { return x*x; }

int kParam;

int calculate_residual()
{
	tt++;
	double eFlow = 0;
      double eVoltage = 0;
	rep(i,0,n-1) b1[i]=0;
      x[0] = 0;
	rep(i,0,n-2)
	{
		x[i + 1] = x[i] + d[i+1].rPath * d[i+1].fPath;
		eFlow += d[i+1].rPath * sqr(d[i+1].fPath);
		eVoltage += sqr(x[i + 1] - x[i]) / d[i+1].rPath;
		b1[i]-=d[i+1].fPath;
		b1[i+1]+=d[i+1].fPath;
      }
	rep(i,0,m-1)
	{
		eFlow += resistance[i+1] * sqr(f[i+1]);
		eVoltage += sqr(x[e[i+1][0]] - x[e[i+1][1]]) / resistance[i+1];
		double fx=(x[e[i+1][0]] - x[e[i+1][1]]) / resistance[i+1];
		b1[e[i+1][0]]+=fx;
		b1[e[i+1][1]]-=fx;
      }
	double xb = 0;
	rep(i,0,n-1) xb += x[i] * b[i];

	double gap = eFlow - (2 * xb - eVoltage);
	fprintf(stderr, "iteration %d: energy = %lf = %lf, xb = %lf, voltage energy = %lf, gap = %.16lf\n", tt, energy, eFlow, xb, eVoltage, gap);
	//mexWarnMsgIdAndTxt("MyToolbox:recursiveLineMex:lol", msg);
	
	double b1norm = 0;
	rep(i,0,n-1) b1norm+=sqr(b1[i]-b[i]);
	b1norm=sqrt(b1norm);
	
	//printf("%.16lf\n",b1norm/bnorm);
	return (b1norm <= double(1e-5) * bnorm) || (LL(tt) > 100000);
}

double scale;

double rS[maxn], sw[maxm];

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
	
	m=0; sw[0]=0;
	rept(it,g.o)
	{
		m++; 
		e[m][0]=get<0>(*it)-1; e[m][1]=get<1>(*it)-1; 
		if (e[m][0]>e[m][1]) swap(e[m][0],e[m][1]);
		resistance[m]=get<2>(*it);
		sw[m]=sw[m-1]+(rS[e[m][1]+1]-rS[e[m][0]+1])/resistance[m]+1;
	}
	//printf("total stretch = %.16lf\n",sw[m]);
	Vec rhs=IO::readMMVec(rhsFile);
	assert(rhs.n==n);
	
	rep(i,0,n-1) b[i]=rhs[i];
	
	bnorm = 0;
	rep(i,0,n-1) bnorm+=sqr(b[i]);
	bnorm = sqrt(bnorm);
	
	double s = 0;
	rep(i,0,n-1)
	{
		s += b[i];
	}
	s = s / double(n);
	rep(i,0,n-1) b[i] -= s;
	
	energy = 0;
	s = 0;
	rep(i,0,n-2)
	{
		s += b[i];
		d[i+1].fPath = -s;
		energy += d[i+1].rPath * sqr(d[i+1].fPath);
	}
	rep(i,1,m) f[i]=0;

	double tsamplecost=0;
	order_num = n*4;
	clock_t t_start = clock();
	while (1)
	{
		clock_t t_start2 = clock();
		rep(i,1,order_num) order[i]=sampleToggle();
		clock_t t_end2 = clock();
		tsamplecost+=double(t_end2 - t_start2) / double(CLOCKS_PER_SEC);
		
		rep(i,1,n-1) d[i].ts=0;
		DS::init(d,n-1);
		
		osum_r[0]=0; osum_rf[0]=0;
		rep(i,1,n-1) 
		{
			osum_r[i]=osum_r[i-1]+d[i].rPath;
			osum_rf[i]=osum_rf[i-1]+d[i].rPath*d[i].fPath;
		}
		
		rep(i,1,order_num) toggle(order[i]);

		DS::dump(ans,n);
		rep(i,1,n-1) d[i].fPath+=ans[i];
		
		if (calculate_residual()) break;
	}
	clock_t t_end = clock();
	double tcost=double(t_end - t_start) / double(CLOCKS_PER_SEC);
	
	Vec vs(n);
	rep(i,0,n-1) vs[i]=x[i];
	Mat A=IO::constructMatrixFromGraph(g);
	double relres=(A*vs-rhs).norm()/rhs.norm();
	//assert(relres<(1e-5)+(1e-9));
	printf("linesolver %0.6e %lld %.3lf %.3lf\n",relres,LL(tt)*order_num,tcost,tsamplecost);
	
	return 0;
}
