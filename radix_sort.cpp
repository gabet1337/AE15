#include <iostream>
#include <vector>
#include <bitset>
#include <algorithm>
#include <time.h>
#include <list>
#include <chrono>
#include <omp.h>
#include <cmath>

#define BASE_BITS 8
#define BASE (1 << BASE_BITS)
#define MASK (BASE-1)
#define DIGITS(v, shift) (((v) >> shift) & MASK)

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

void omp_lsd_radix_sort(size_t n, unsigned data[]) {
  unsigned * buffer = (unsigned*)malloc(n*sizeof(unsigned));
  int total_digits = sizeof(unsigned)*8;
 
  //Each thread use local_bucket to move data
  size_t i;
  for(int shift = 0; shift < total_digits; shift+=BASE_BITS) {
    size_t bucket[BASE] = {0};
 
    size_t local_bucket[BASE] = {0}; // size needed in each bucket/thread
    //1st pass, scan whole and check the count
#pragma omp parallel firstprivate(local_bucket)
    {
#pragma omp for schedule(static) nowait
      for(i = 0; i < n; i++){
        local_bucket[DIGITS(data[i], shift)]++;
      }
#pragma omp critical
      for(i = 0; i < BASE; i++) {
        bucket[i] += local_bucket[i];
      }
#pragma omp barrier
#pragma omp single
      for (i = 1; i < BASE; i++) {
        bucket[i] += bucket[i - 1];
      }
      int nthreads = omp_get_num_threads();
      int tid = omp_get_thread_num();
      for(int cur_t = nthreads - 1; cur_t >= 0; cur_t--) {
        if(cur_t == tid) {
          for(i = 0; i < BASE; i++) {
            bucket[i] -= local_bucket[i];
            local_bucket[i] = bucket[i];
          }
        } else { //just do barrier
#pragma omp barrier
        }
 
      }
#pragma omp for schedule(static)
      for(i = 0; i < n; i++) { //note here the end condition
        buffer[local_bucket[DIGITS(data[i], shift)]++] = data[i];
      }
    }
    //now move data
    unsigned* tmp = data;
    data = buffer;
    buffer = tmp;
  }
  free(buffer);
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
    cout << A[i] << " goes to bucket3 " << d << " at position " << next3[d] << endl;
    next3[d]++;
  }
  //omp barrier
  //#pragma omp barrier 
  //calculate prefix sum
  vui output_indices(256,0);
  output_indices[0] = 0;
  for (int i = 1; i < 256; i++) output_indices[i] = next3[i] + output_indices[i-1];
  // for (int i = 0; i < 256; i++) cout << output_indices[i] << " "; 
  // cout << endl;
  //sort each bucket in parallel using counting sort from the LSD
  //for each bucket
  //#pragma omp parallel for
  for (size_t i = 0; i < 256; i++) {
    size_t bucket_i_size = next3[i];
    vvui bucket0(256, vui());
    for (size_t j = 0; j < 256; j++) bucket0[j].reserve(bucket_i_size);
    vui next0(256,0);
    //for each item in bucket i on first digit
    for (size_t j = 0; j < bucket_i_size; j++) {
      unsigned int d = (bucket3[i][j] & 0xFF);
      cout << bucket3[i][j] << " goes to bucket0 " << d << " at position " << next0[d] << endl;
      bucket0[d][next0[d]] = bucket3[i][j];
      next0[d]++;
    }

    vvui bucket1(256,vui());
    for (size_t j = 0; j < 256; j++) bucket1[j].reserve(bucket_i_size);
    vui next1(256,0), histogram2(256,0), histogram1(256,0);
    for (size_t j = 0; j < 256; j++) {
      // bucket0 = bucket0[j]
      size_t bucket0_j_size = next0[j];
      for (size_t k = 0; k < bucket0_j_size; k++) {
        unsigned int d = ((bucket0[j][k] & 0xFF00) >> 8);
        bucket1[d][next1[d]] = bucket0[j][k];
        cout << bucket0[j][k] << " goes to bucket1 " << d << " at position  " << next1[d] << endl;
        next1[d]++;
        histogram1[d]++;
        d = ((bucket0[j][k] & 0xFF0000) >> 16);
        histogram2[d]++;
      }
    }

    //for each bucket in bucket1
    for (int b = 255; b >= 0; b--) {
      //for each item in bucket1
      size_t bucket1_b_size = next1[b];
      for (size_t j = 0; j < bucket1_b_size; j++) {
        unsigned int d0 = (bucket1[b][j] & 0xFF);
        unsigned int d1 = (bucket1[b][j] & 0xFF00) >> 8;
        unsigned int d2 = (bucket1[b][j] & 0xFF0000) >> 16;
        unsigned int d3 = (bucket1[b][j] >> 24);
        cout << bucket1[b][j] << " d: " << d3 << " " << d2 << " " << d1 << " " << d0 << endl;
        cout << output_indices[d3] << " " << histogram2[d2] << " " << histogram1[d1] << " " << next0[d0] << endl;
        int k = output_indices[d3] + (histogram2[d2] - histogram1[d2]);
        cout << k << " = " << output_indices[d3] << " - " << histogram2[d2] << endl;
        A[k] = bucket1[b][j];
        // cout << k << " " << bucket1[b][j] << endl; 
      }
    }
  }
}

int main() {
  srand(time(NULL));
  size_t size = 100000000*2;
  unsigned *data = (unsigned*)malloc(size * sizeof(unsigned));
  vector<unsigned int> data2(size,0);
  for (size_t i = 0; i < size; i++) {
    int r = rand();
    data[i] = r;
    data2[i] = r;
  }


  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  //parallel_radix_sort(s);
  omp_lsd_radix_sort(size, data);
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
  cout << duration << endl;  

  t1 = high_resolution_clock::now();
  radix_sort(data2, 8);  
  t2 = high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
  cout << duration << endl;  


  // cout << size << endl;
  // for (size_t i = 0; i < size; i++)
  //   cout << data[i] << " ";
  // cout << endl;
  if (is_sorted(data, data+size)) cout << "SUCCESS!!!" << endl;
  else cout << "FAILED!!!" << endl;
  if (is_sorted(data2.begin(), data2.end())) cout << "SUCCESS!!!" << endl;
  else cout << "FAILED!!!" << endl;
    

  

  return 0;
}
