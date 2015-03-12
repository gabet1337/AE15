#include "counting_sort.cpp"
#include <iostream>
#include <vector>
#include <bitset>
#include <algorithm>
#include <time.h>
#include <list>
#include <chrono>
#include <omp.h>
#include <cmath>
using namespace std;
using namespace std::chrono;

typedef vector<unsigned int> vui;
typedef vector<vui> vvui;

//k should divide 32 but should not be 32. (1,2,4,8,16)
void radix_sort(vector<unsigned int> &A, unsigned int k) {
  unsigned int mask = (1 << k) - 1;
  size_t mask_size = mask+1;
  for (unsigned int i = 0; i < 32/k; i++) {
    mask = mask << (k * (i != 0));
    //counting sort using the mask
    vector<unsigned int> B, C(mask_size, 0);
    B.reserve(A.size());
    for (size_t x = 0; x < A.size(); x++)
      C[(A[x]&mask)>>(k*i)]++;
    for (int j = 1; j < (int)mask_size; j++)
      C[j] += C[j-1];
    for (int j = (int)A.size()-1; j >= 0; j--)
      B[--C[(A[j]&mask)>>(k*i)]] = A[j];
    A = B;
  }
}

void parallel_radix_sort(vui &A) {
  size_t size = A.size();
  vui output;
  output.reserve(size);
  vvui bucket3(256, vui()); //TODO: look at initializing in parallel
  for (int i = 0; i < 256; i++) bucket3[i].reserve(size);
  vui next3(256,0);
  //sort the digits by MSD using counting sort in parallel
  //#pragma omp parallel for
  for (size_t i = 0; i < size; i++) {
    unsigned int d = (A[i] >> 24);
    bucket3[d][next3[d]] = A[i];
    next3[d]++;
    cout << d << " -> " << A[i] << " = " << next3[d] << endl;
  }
  //omp barrier
  #pragma omp barrier
  //calculate prefix sum
  vui output_indices(256,0);
  output_indices[0] = next3[0]-1;
  for (int i = 1; i < 256; i++) output_indices[i] = next3[i] + output_indices[i-1];
  for (int i = 0; i < 256; i++) cout << output_indices[i] << " ";
  cout << endl;
  //sort each bucket in parallel using counting sort from the LSD
  //for each bucket
  //#pragma omp parallel for
  for (size_t i = 0; i < 256; i++) {
    vvui bucket0(256, vui());
    for (size_t j = 0; j < 256; j++) bucket0[j].reserve(next3[i]);
    vui next0(256,0);
    //for each item in bucket i on first digit
    for (size_t j = 0; j < next3[i]; j++) {
      unsigned int d = (bucket3[i][j] & 0xFF);
      bucket0[d][next0[d]] = bucket3[i][j];
      next0[d]++;
    }

    vvui bucket1(256,vui());
    for (size_t j = 0; j < 256; j++) bucket1[j].reserve(next3[i]);
    vui next1(256,0), histogram(256,0);
    for (size_t j = 0; j < 256; j++) {
      // bucket0 = bucket0[j]
      for (size_t k = 0; k < next0[j]; k++) {
        unsigned int d = ((bucket0[j][k] & 0xFF00) >> 8);
        bucket1[d][next1[d]] = bucket0[j][k];
        next1[d]++;
        d = ((bucket0[j][k] & 0xFF0000) >> 16);
        histogram[d]++;
        cout << d << " --> " << bucket0[j][k] << " = " << histogram[d] << endl; 
      }
    }

    //for each bucket in bucket1
    for (size_t b = 0; b < 256; b++) {
      //for each item in bucket1
      for (size_t j = 0; j < next1[b]; j++) {
        unsigned int d2 = (bucket1[b][j] & 0xFF0000) >> 16;
        unsigned int d3 = (bucket1[b][j] >> 24);
        int k = output_indices[d3] - histogram[d2];
        histogram[d2]--;
        A[k] = bucket1[b][j];
        cout << k << " " << bucket1[b][j] << endl;
      }
    }
  }
}

int main() {
  srand(time(NULL));

  vector<unsigned int> s;
  s.push_back(0);
  s.push_back(0xFF);
  s.push_back(0xFF00);
  s.push_back(0xFF0000);
  s.push_back(0xFF000000);

  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  parallel_radix_sort(s);
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
  cout << duration << endl;  
  
  cout << s.size() << endl;
  for (size_t i = 0; i < s.size(); i++)
    cout << s[i] << " ";
  cout << endl;  
  if (is_sorted(s.begin(), s.end())) cout << "SUCCESS!!!" << endl;
  else cout << "FAILED!!!" << endl;
  

  return 0;
}
