////////data generator
//   data format (note: 0 indexing):
//     line 1: n (#V), m (#E), k (# cycle toggles)
//     lines 2 - n: tree edges, in u v w format
//     lines n + 1 ... m + 1: off tree edges, in u v r format:
//       vert1, vert2, weight = 1 / r
//     lines m + 2 ... 2m + 1: initial flow on edges
//     lines 2m + 2... 2m + 1 + k: one id between n and m each, the id of the cycles toggled

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

int n, m, seed;

const int MAXN = 1<<24;
const int RNG = 1000;
double diag[1000000];

D r[MAXN];
D rS[MAXN];

set<pair<int,int> > has;

int main(int argt, char **args) {

  sscanf(args[1], "%d", &n);
  sscanf(args[2], "%d", &m);

  printf("%%%%MatrixMarket matrix coordinate real symmetric\n%%\n");
  cout << n << ' ' << n  << ' ' << m+n+n << endl;

  rS[0] = 0;
  FOR(i, n - 1) {
    r[i] = rand() % 1000+1;
    printf("%d %d %0.6e\n", i+1, i + 2, -1./r[i]);
    diag[i]+=1./r[i];
    diag[i+1]+=1./r[i];
    rS[i + 1] = rS[i] + r[i];
  }
  has.clear();
  FOR(i, m - n + 1) {
    int u = rand() % n;
    int v = rand() % n;
    if(u > v) swap(u, v);
    while(u == v || v == u + 1 || has.find(MP(u, v)) != has.end()) {
      u = rand() % n;
      v = rand() % n;
      if(u > v) swap(u, v);
    }
    has.insert(MP(u, v));   
    printf("%d %d %.4lf\n", u+1, v+1, -1./(rS[v]-rS[u]));
    diag[u]+=1./(rS[v]-rS[u]);
    diag[v]+=1./(rS[v]-rS[u]);
  }


  for(int i=0; i < n; ++i) {
    printf("%d %d %0.6e\n",i+1,i+1,diag[i]);
  }


  return 0;
}

