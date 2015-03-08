#include <algorithm>
#include <vector>
#include <iostream>
using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int,int> ii;
typedef pair<ii, ii> dim;
//A: m x n
//B: n x p
//C: m x p
vvi rec_mult(vvi A, vvi B, vvi C) {
  int m = A.size();
  int n = B.size();
  int p = A[0].size();

  if (m == 1 && n == 1 && p == 1) {
    return C + A[0][0] * B[0][0];
  } else if (m > max(n,p)) {
    vvi C1(C.size()/2, vi(C[0].size(), 0));
    for (int i = 0; i < C.size()/2; i++)
      for (int j = 0; j < C[0].size(); j++)
        C1[i][j] = C[i][j];
    vvi C2(C.size()/2+(C.size()%2), vi(C[0].size(), 0));
    for (int i = C.size()/2+(C.size()%2); i < C.size(); i++)
      for (int j = 0; j < C[0].size(); j++)
        C1[i][j] = C[i][j];

    

  } else if (n > max(m,p)) {

  } else if (p > max(m,n)) {
    
  }
}

dim get_dim(vvi &x) {
  return dim(ii(0,0), ii(x.size(), x[0].size()));
}

int main() {

  vvi A,B,C;
  A.assign(3, vi(3,0));
  B.assign(3, vi(3,0));
  C.assign(3, vi(3,0));
  rec_mult(A,B,C,get_dim(A),get_dim(B),get_dim(C));
  for (int i = 0; i < C.size(); i++) {
    for (int j = 0; j < C[0].size(); j++)
      cout << C[i][j] << " ";
    cout << endl;
  }
  return 0;
}
