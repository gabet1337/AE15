#include "counting_sort.cpp"
#include <iostream>
#include <vector>
#include <bitset>
#include <algorithm>
#include <time.h>
#include <chrono>
using namespace std;
using namespace std::chrono;
//k should divide 32 but should not be 32. (1,2,4,8,16)
void radix_sort(vector<unsigned int> &A, int k) {
  unsigned int mask = (1 << k) - 1;
  size_t mask_size = mask+1;
  for (unsigned int i = 0; i < 32/k; i++) {
    mask = mask << (k * (i != 0));
    //counting sort using the mask
    vector<unsigned int> B((int)A.size(), 0), C(mask_size, 0);
    for (size_t j = 0; j < A.size(); j++)
      ++C[(A[j]&mask)>>(k*i)];
    for (int j = 1; j < mask_size; j++)
      C[j] += C[j-1];
    for (int j = (int)A.size()-1; j >= 0; j--)
      B[--C[(A[j]&mask)>>(k*i)]] = A[j];
    A = B;
  }
}

int main() {
  srand(time(NULL));
  int k = 4;

  vector<unsigned int> s;
  for (int i = 0; i < 10000000; i++) {
    int r = rand(); 
    s.push_back(r);
  }
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  radix_sort(s, 1);
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
  cout << duration << endl;  
  t1 = high_resolution_clock::now();
  radix_sort(s, 2);
  t2 = high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
  cout << duration << endl;
  t1 = high_resolution_clock::now();
  radix_sort(s, 4);
  t2 = high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
  cout << duration << endl;
  t1 = high_resolution_clock::now();
  radix_sort(s, 8);
  t2 = high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
  cout << duration << endl;
  t1 = high_resolution_clock::now();
  radix_sort(s, 16);
  t2 = high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
  cout << duration << endl;

  


  // for (int i = 0; i < s.size(); i++)
  //   cout << s[i] << " ";
  // cout << endl;
  if (is_sorted(s.begin(), s.end())) cout << "SUCCESS!!!" << endl;
  else cout << "FAILED!!!" << endl;
  

  return 0;
}
