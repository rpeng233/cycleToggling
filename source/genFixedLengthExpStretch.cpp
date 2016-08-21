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

int n, m, hop, k, seed;

const int MAXN = 1<<24;
const int RNG = 1000;

D r[MAXN];
D rS[MAXN];

set<pair<int,int> > has;

int main(int argt, char **args) {
  fprintf(stderr, "USAGE: #vertices on path, #hopCount, #toggles, randomness seed\n");
  fflush(stdout);
// each path connects vertices hop apart
// but still have stretch 1
  sscanf(args[1], "%d", &n);
  sscanf(args[2], "%d", &hop);
  if(hop > n) {
    fprintf(stderr, "hop length of paths need to be < n");
    return 0;
  }

  sscanf(args[3], "%d", &k);
  sscanf(args[4], "%d", &seed);
  srand(seed);
  m = n - 1 + (n - hop);
  printf("%d %d %d\n", n, m, k);
  rS[0] = 0;
  FOR(i, n - 1) {
    r[i] = rand() % 10000 + 1;
    printf("%d %d %.16lf\n", i, i + 1, r[i]);
    rS[i + 1] = rS[i] + r[i];
  }

  has.clear();
  FOR(i, n - hop) {
    int u = i;
    int v = i + hop; 
    has.insert(MP(u, v));   
    double r=rS[v]-rS[u];
    while (double(rand())/double(RAND_MAX)<0.5) r/=2;
    printf("%d %d %.16lf\n", u, v, r);
  }
  //FOR(i, n) {
  //  printf("%d.0\n", rand() % RNG);
  //}
  FOR(i,m) {
    if(i==0)
      printf("%d.0\n", -1);
    else if(i==(n-1))
      printf("%d.0\n", 1);
    else
      printf("%d.0\n", 0);
      }
  /*   FOR(i, k) {
    printf("%d\n", rand() % (m - n + 1) + n - 1);
    }*/
  return 0;
}

