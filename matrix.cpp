#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

const int N=10000;

int M[N][N];

int main() {
  clock_t start,end;
  
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      M[i][j]=i+j+1;

  for (int r=0; r<5; r++) {
    start = clock();
    int sum=0;
    for (int j=0; j<N; j++)
      for (int i=0; i<N; i++)
	sum+=M[i][j];
    
    end = clock();
    double time = (end - start) / (double)(CLOCKS_PER_SEC / 1000);
    cout << "column order: " << time/1000 << " seconds" << "  sum=" << sum << endl;
    
    start = clock();
    
    sum=0;
    for (int i=0; i<N; i++)
      for (int j=0; j<N; j++)
	sum+=M[i][j];
    
    end = clock();
    time = (end - start) / (double)(CLOCKS_PER_SEC / 1000);
    cout << "row order: " << time/1000 << " seconds"  << "  sum=" << sum << endl;
  };    
};
