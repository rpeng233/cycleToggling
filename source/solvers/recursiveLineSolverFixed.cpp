#include <cassert>
#include <cstdio> 
#include <cstring> 
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <algorithm> 
#include <string> 
#include <vector> 
#include <set>
#include <map> 
#include <queue>
#include <complex>
#include "../util/io.h"
#include "../util/graph.h"
using namespace std; 
typedef long long ll;
typedef double D;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int,int> ii;
typedef vector<ii> vii;
typedef vector<vii> vvii;

#define MP make_pair 
#define A first 
#define B second 

#define PB push_back 
#define FR(i, a, b) for(int i=(a); i<(b); i++) 
#define FOR(i, n) FR(i, 0, n) 
#define RF(i, a, b) for(int i=(b)-1; i>=(a); i--) 
#define ROF(i, n) RF(i, 0, n) 
#define EACH(it,X) for(__typeof((X).begin()) it=(X).begin(); it!=(X).end(); ++it) 


const int MAXN = 5000000;
const int MAXM = 10000010;
const int CHUNK= 4;

D energy;
D b[MAXN], x[MAXN];
D b1[MAXN];

D sqr(D x) {
  return x * x;
}

int n, m, k;
D path[MAXN * 2][3];

struct updateType{
  int u, v, id;
  bool operator < (const updateType &o) const {
    return MP(u, v) < MP(o.u, o.v);
  }
} update[MAXN * 2];

struct edgeType {
  int u, v;
  D r, f;
  bool operator < (const edgeType &o) const {
    return MP(u, v) < MP(o.u, o.v);
  }
} e[MAXN];

const int WORD_LENGTH = 32;
unsigned int mask[MAXN * 2];
int maskSum[MAXN];
int offset; //offset into bitmask array

//BITMASK, with SUMMATION
inline void setBit(int id) {
  mask[offset + id / WORD_LENGTH] |= (1 << (id % WORD_LENGTH));
}

inline int getBit(int id) {
  return (mask[offset + id / WORD_LENGTH] >> (id % WORD_LENGTH)) % 2;
}

inline void buildSum(int len) {
  int s = 0;
  FOR(i, len) {
    maskSum[i] = s;
    s += __builtin_popcount(mask[offset + i]); 
  }
}

inline int getSumBit(int id) {
  return maskSum[id / WORD_LENGTH] +
          __builtin_popcount(((1 << id) - 1) & mask[offset + id / WORD_LENGTH]);
}

int restriction(int pathL, int pathR, int updateL1, int updateR1, int updateR) {
//clock_t tRestriction_start = clock();

//clock_t tRemove_start = clock();
  int n = pathR - pathL + 1;
  int nSmall = n / WORD_LENGTH + 1;
  FOR(i, nSmall) {
    mask[offset + i] = 0;
  }
  FR(i, updateL1, updateR1) {
    setBit(update[i].u);
    setBit(update[i].v);
  }
//clock_t tRemove_end = clock();
//if(n > 2000000) printf("time for figuring out what to remove: %0.3lfs\n", D(tRemove_end - tRemove_start) / D(CLOCKS_PER_SEC));
//clock_t tPathCompress_start = clock();
  int n1 = 0;
  D s0 = 0;
  D s1 = 0;
  int i = 0;


  FOR(i1, nSmall) {
    unsigned int m1 = mask[offset + i1];
    FOR(i2, WORD_LENGTH) {
//printf("keeping %d\n", i);
      if(i == n) {
        break;
      }
//printf("%d    %d %d   %d %d    %d\n", i, i1, i2, getBit(i), m1, n);
//assert(getBit(i) == (m1 >> i2) % 2);
      if(m1 % 2 == 1) {//if(getBit(i)) {
        if(n1 > 0) {
          path[pathR + n1 - 1][0] = s0;
          path[pathR + n1 - 1][1] = s1;
          path[pathR + n1 - 1][2] = 0;
        }
        n1++;
        s0 = 0;
        s1 = 0;
      }
      s0 += path[pathL + i][0];
      s1 += path[pathL + i][1];
      i++;
      m1 = m1 >> 1;
    }
  }

//clock_t tPathCompress_end = clock();
//if(n > 2000000) printf("time for pathCompress: %0.3lfs\n", D(tPathCompress_end - tPathCompress_start) / D(CLOCKS_PER_SEC));

//clock_t tRelabel_start = clock();
  buildSum(nSmall);
  FR(i, updateL1, updateR1) {
//printf("...%d, %d\n", updateR + i - updateL1, update[i].id);
    update[updateR + i - updateL1].u = getSumBit(update[i].u);
    update[updateR + i - updateL1].v = getSumBit(update[i].v);
//taking out memory access for profiling purposes
//update[updateR + i - updateL1].u = rand() % n1;
//update[updateR + i - updateL1].v = rand() % n1;
    update[updateR + i - updateL1].id = update[i].id;
  }
//clock_t tRelabel_end = clock();
//if(n > 2000000) printf("time for relabeling: %0.3lfs\n", D(tRelabel_end - tRelabel_start) / D(CLOCKS_PER_SEC));
 
//clock_t tRestriction_end = clock();
//if(n > 2000000) {printf("time for restriction: %0.3lfs\n", D(tRestriction_end - tRestriction_start) / D(CLOCKS_PER_SEC)); printf("end of restriction, %d->%d\n", n, n1); }
  return n1;
}

void prolongation(int pathL, int pathR, int newLength) {
  int n = pathR - pathL + 1;
  int nSmall = n / WORD_LENGTH + 1;

  int pre;
  int num = 0;
  int i = 0;
  FOR(i1, nSmall) {
    unsigned int m1 = mask[offset + i1];
    FOR(j, WORD_LENGTH) {
      if(i == n) {
        break;
      }
      if(m1 % 2 == 1) {
//printf("prolongating from: %d\n", i1);
        num++;
        pre = i;
      }
      if(num > 0 && num < newLength) {
        path[pathL + i][1] += path[pathL + i][0] * path[pathR + num - 1][2];
        path[pathL + i][2] += path[pathR + num - 1][2];
      }
      i++;
      m1 = m1 >> 1;
    }
  }
}

void solve(int pathL, int pathR, int updateL, int updateR) {
//printf("offset = %d, path[%d,%d], update[%d,%d]\n", offset, pathL, pathR, updateL, updateR);
//FR(i, pathL, pathR) printf("r = %lf, rf = %lf, shift = %lf\n", path[i][0], path[i][1], path[i][2]);
//FR(i, updateL, updateR) printf("%d %d %d\n", update[i][0], update[i][1], update[i][2]);

// switching to brute force since 50n is ok: should be || here
// if(updateR - updateL <= (1<<6)) return;
  if(updateR - updateL <= 50 || pathR - pathL <= 50) {
    FR(i, updateL, updateR) {
      int u = update[i].u;
      int v = update[i].v;
      int id = update[i].id;
//printf("%d:: %d %d\n", id, u, v);
      D sum_rf = -e[id].f * e[id].r;
      D sum_r = e[id].r;
      FR(j, u, v) {
        sum_r += path[pathL + j][0];
        sum_rf += path[pathL + j][1];
      }
      D delta = sum_rf / sum_r;
//printf("%lf %lf, amount moved = %lf\n", sum_rf, sum_r, delta);
      e[id].f += delta;
      FR(j, u, v) {
        path[pathL + j][1] -= delta * path[pathL + j][0];
        path[pathL + j][2] -= delta;
      }
      energy -= delta * sum_rf;
    }
    return;
  }
  int n = pathR - pathL + 1;
  int nSmall = n / WORD_LENGTH + 1;

  int chunkSize = (updateR - updateL) / CHUNK;
  FOR(i, CHUNK) {
    int updateL1 = updateL + i * chunkSize;
    int updateR1 = updateL1 + chunkSize;
    if(i == CHUNK - 1) {
      updateR1 = updateR;
    }
    int n1 = restriction(pathL, pathR, updateL1, updateR1, updateR);
    offset += nSmall;
    solve(pathR, pathR + n1 - 1, updateR, updateR + (updateR1 - updateL1));
    offset -= nSmall;
    prolongation(pathL, pathR, n1);
//puts("===result from prolongation");
//FR(i, pathL, pathR) printf("r = %lf, rf = %lf\n", path[i][0], path[i][1]);
  }
}

double rS[MAXN], sw[MAXM];

int sampleToggle()
{
	double x=double(rand())/double(RAND_MAX)*sw[m];
	int t=lower_bound(sw+1,sw+m+1,x)-sw;
	return t-1;
}


int main(int argc, char *argv[])
{
	if (argc!=3) { printf("Usage: ./linesolver graphFile rhsFile\n"); return 0; }
	
	string graphFile=argv[1];
	string rhsFile=argv[2];
	
	GraphSP g=IO::readGraphSP(graphFile);
	n=g.n;
  ios::sync_with_stdio(0);
  int u, v;
  D w;
  rS[1]=0;
  rep(i,1,n-1)
  {
	  assert(g.e[i].size()==1);
	assert(g.e[i][0].first==i+1);
    path[i-1][0] =g.e[i][0].second;
    rS[i+1]=rS[i]+path[i-1][0];  
  }
  m=g.o.size();
  sw[0]=0;
  FOR(i, m) {
    int i1 = i;
    e[i].u=get<0>(g.o[i1]) - 1;
    e[i].v=get<1>(g.o[i1]) - 1;
    e[i].r=get<2>(g.o[i1]);

//if(e[i].u == 2 || e[i].v == 2) printf("!@#$!@#$ %d %d\n", e[i].u, e[i].v);
    if(e[i].u > e[i].v) {
      swap(e[i].u, e[i].v);
    }
    sw[i+1]=sw[i]+(rS[e[i].v+1]-rS[e[i].u+1])/e[i].r+1;
  }
  Vec rhs=IO::readMMVec(rhsFile);
assert(rhs.n==n);
rep(i,0,n-1) b[i]=rhs[i];
  double s = 0;
  FOR(i, n) {
    s += b[i];
  }
  s = s / double(n);
  FOR(i, n) {
    b[i] -= s;
  }
	
  clock_t t_start = clock();
  energy = 0;
  s = 0;
  FOR(i, n - 1) {
    s += b[i];
    path[i][1] = -s * path[i][0];
    path[i][2] = 0;
//printf("%g %g\n", sqr(fPath[i]) * rPath[i], rPath[i]);
    energy += path[i][0] * sqr(s);
  }
  FOR(i, m) {
    e[i].f = 0;
  }
//fprintf(stderr, ":::initial energy = %lf\n", energy);

  D bNorm = 0;
  FOR(i, n) {
    bNorm += sqr(b[i]);
  }
  bNorm = sqrt(bNorm);


  Vec ans(n);
  Mat A=IO::constructMatrixFromGraph(g);

  rep(i,0,n-1) {
    ans[i] = b[i];
  }
//  printf("initial norms: %lf %lf\n", bNorm, ans.norm());

  double tsamplecost=0;
  int order_num=n;
int tt;
  for(tt = 0;;++tt) {
    int k = order_num;
    clock_t t_start2 = clock();
    FOR(i, k) {
      update[i].id = sampleToggle();
      update[i].u = e[update[i].id].u;
      update[i].v = e[update[i].id].v;
    }
    clock_t t_end2 = clock();
    tsamplecost+=double(t_end2 - t_start2) / double(CLOCKS_PER_SEC);
    offset = 0;
    solve(0, n - 1, 0, k);

    double eFlow = 0;
    double eVoltage = 0;
    x[0] = 0;
    FOR(i, n) {
      b1[i] = 0;
    }

    FOR(i, n - 1) {
      x[i + 1] = x[i] + path[i][1];
      eFlow += path[i][0] * sqr(path[i][1] / path[i][0]);
      eVoltage += sqr(x[i + 1] - x[i]) / path[i][0];
      D f1 = path[i][1] / path[i][0];
      b1[i] -= f1;
      b1[i + 1] += f1;
    }

//printf("========only path=====\n%0.10lf\n", b1[1]);
//printf("%0.10lf\n", (x[1] - x[0]) / D(394)
//     + (x[1] - x[2]) / D(14));
 
    FOR(i, m) {
      eFlow += e[i].r * sqr(e[i].f);
      eVoltage += sqr(x[e[i].u] - x[e[i].v]) / e[i].r;
      D f1 = (x[e[i].u] - x[e[i].v]) / e[i].r;
//if(e[i].u == 1 || e[i].v == 1) printf("!@#$!@#$ %d %d", e[i].u, e[i].v);
      b1[e[i].u] += f1;
      b1[e[i].v] -= f1;
    }
    double xb = 0;
    FOR(i, n) {
      xb += x[i] * b[i];
    }
    double gap = eFlow - (2 * xb - eVoltage);

    D b1Norm = 0;
    FOR(i, n) {
      b1Norm += sqr(b1[i] - b[i]);
//printf("%lf %lf\n", b[i], b1[i]);
    }
    b1Norm = sqrt(b1Norm);

    rep(i,0,n-1) {
      ans[i] = x[i];
    }
    ans = A * ans;
/*
FOR(i, 10) printf("WTF: %g  %d\n", fabs(b1[i] - ans[i]), g.e[i].size());
printf("===%d\n", A.values.size());
printf("%d      mine: %0.10lf matrix's: %0.10lf\n", g.e[1].size(), b1[1], ans[1]);
//FOR(j, g.e[1].size())  printf("edge: to %d, resistance %0.10lf\n", (g.e[1])[j].first, (g.e[1])[j].second);
EACH(ii, A.values) if((ii -> x == 1 || ii -> y == 1) && (ii -> x != ii -> y)) printf("%d %d %0.10lf\n", ii -> x, ii -> y, D(1) / (-ii -> z));
D v1 = (x[1] - x[0]) / D(394)
     + (x[1] - x[2]) / D(14)
     + (x[1] - x[141275]) / D(70423716)
     + (x[1] - x[457061]) / D(228338004);
printf("actual::: %0.10lf\n", v1);
*/
//printf("?????%0.10lf %0.10lf\n", b1Norm, (ans - rhs).norm());
//exit(0); 
    fprintf(stderr, "iteration %d: residue = %lf, energy = %lf = %lf, xb = %lf, voltage energy = %lf, gap = %lf\n", tt, b1Norm / bNorm, energy, eFlow, xb, eVoltage, gap);
    //printf("residual %0.6e\n",b1Norm/bNorm);
    if(b1Norm / bNorm < double(1e-5)) {
      break;
    }
    if(LL(tt) > 100000) {
      break;
    }
  }

  rep(i,0,n-1) {
    ans[i] = x[i];
  }
  double relres=(A*ans-rhs).norm()/rhs.norm();
  //  printf("%.16lf\n",relres);
  //assert(relres<(1e-5)+(1e-9));

  clock_t t_end = clock();
  double tcost=double(t_end - t_start) / double(CLOCKS_PER_SEC);
  printf("reclinesolver %0.6e %lld %.3lf %.3lf\n",relres,LL(tt)*order_num,tcost,tsamplecost);
  
/*
  puts("Outputting resulting vector to x.txt");
  FILE *f_out = fopen("x.txt", "w");
  FOR(i, n) {
    fprintf(f_out, "%0.3g\n", x[i]);
  }
  fclose(f_out);
  */
}
