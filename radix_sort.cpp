#include <iostream>
#include <vector>
#include <bitset>
#include <algorithm>
#include <time.h>
#include <list>
#include <chrono>
#include <omp.h>
#include <cmath>
#include <fstream>
#include <random>
#include <limits.h>
#include <papi.h>
#define BASE_BITS 8
#define BASE (1 << BASE_BITS)
#define MASK (BASE-1)
#define DIGITS(v, shift) (((v) >> shift) & MASK)
#define BUFFER_BASE 32

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

void test_running_time() {
  size_t num_experiments = 10;
  long long single_core_average = 0, multi_core_average = 0;
  ofstream single, multi;
  single.open("radix_single_running_time.dat");
  multi.open("radix_multi_running_time.dat");
  single << "#x\ty" << endl;
  multi << "#x\ty" << endl;
  for (double s = 20; s <= 28; s+=0.5) {
    size_t size = (size_t)pow((double)2, s);
    for (size_t e = 0; e < num_experiments; e++) {
      cout << "Generating testdata for input size: " << size << endl;
      unsigned *data = (unsigned*)malloc(size * sizeof(unsigned));
      vector<unsigned int> data2(size,0);
      random_device rd;
      mt19937 gen(rd());
      uniform_int_distribution<unsigned int> dis(0, UINT_MAX);
      for (size_t i = 0; i < size; i++) {
        unsigned int r = dis(gen);
        data[i] = r;
        data2[i] = r;
      }
      high_resolution_clock::time_point t1 = high_resolution_clock::now();
      omp_lsd_radix_sort(size,data);
      high_resolution_clock::time_point t2 = high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
      multi_core_average += duration;
      t1 = high_resolution_clock::now();
      radix_sort(data2,BASE_BITS);
      t2 = high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
      single_core_average += duration;
      delete(data);
      data2.clear();
    }
    multi << size << "\t" << multi_core_average/num_experiments << endl;
    single << size << "\t" << single_core_average/num_experiments << endl;
  }
  single.close();
  multi.close();
}
typedef pair<size_t, long long> uill;
typedef vector<uill> vuill;
typedef pair<vuill, vuill> pvuill;
typedef pair<pvuill,pvuill> ppvuill;
ppvuill test_papi_events(int events[]) {
  size_t num_experiments = 10;
  PAPI_library_init(PAPI_VER_CURRENT);
  long long values[2] = {(long long)0};

  long long val11 = 0;
  long long val12 = 0;
  long long val21 = 0;
  long long val22 = 0;
  pvuill exp1,exp2;

  for (double s = 5; s <= 28; s+=0.5) {
    size_t size = (size_t)pow((double)2, s);
    for (size_t e = 0; e < num_experiments; e++) {
      cout << "Generating testdata for input size: " << size << endl;
      unsigned *data = (unsigned*)malloc(size * sizeof(unsigned));
      vector<unsigned int> data2(size,0);
      random_device rd;
      mt19937 gen(rd());
      uniform_int_distribution<unsigned int> dis(0, UINT_MAX);
      for (size_t i = 0; i < size; i++) {
        unsigned int r = dis(gen);
        data[i] = r;
        data2[i] = r;
      }
      cout << "Running experiment..." << endl;
      
      PAPI_start_counters(events, 2);
      radix_sort(data2,BASE_BITS); 
      PAPI_stop_counters(values,2);
      val11 += values[0];
      val12 += values[1];

      PAPI_start_counters(events,2);
      omp_lsd_radix_sort(size,data);
      PAPI_stop_counters(values,2);
      val21 += values[0];
      val22 += values[1];
      
      delete(data);
      data2.clear();
    }

    exp1.first.push_back(make_pair(size, val11/num_experiments));
    exp1.second.push_back(make_pair(size, val12/num_experiments));
    exp2.first.push_back(make_pair(size, val21/num_experiments));
    exp2.second.push_back(make_pair(size, val22/num_experiments));
    
  }
  return make_pair(exp1,exp2);
}

void test_L2() {

  int events[2] = {PAPI_L2_TCM, PAPI_L2_TCA};
  ppvuill res = test_papi_events(events);

  ofstream single,multi;
  single.open("radix_single_L2_ratio.dat");
  multi.open("radix_multi_L2_ratio.dat");
  single << "#x\ty" << endl;
  multi << "#x\ty" << endl;
  vuill single_tcm = res.first.first;
  vuill single_tca = res.first.second;
  vuill multi_tcm = res.second.first;
  vuill multi_tca = res.second.second;
  //single
  for (size_t i = 0; i < single_tcm.size(); i++) {
    single << single_tcm[i].first << "\t" << (double)single_tcm[i].second/(double)single_tca[i].second << endl;
  }
  for (size_t i = 0; i < multi_tcm.size(); i++) {
    multi << multi_tcm[i].first << "\t" << (double)multi_tcm[i].second/(double)multi_tca[i].second << endl;
  }
  single.close();
  multi.close();
}

void test_L3() {

  int events[2] = {PAPI_L3_TCM, PAPI_L3_TCA};
  ppvuill res = test_papi_events(events);

  ofstream single,multi;
  single.open("radix_single_L3_ratio.dat");
  multi.open("radix_multi_L3_ratio.dat");
  single << "#x\ty" << endl;
  multi << "#x\ty" << endl;
  vuill single_tcm = res.first.first;
  vuill single_tca = res.first.second;
  vuill multi_tcm = res.second.first;
  vuill multi_tca = res.second.second;
  //single
  for (size_t i = 0; i < single_tcm.size(); i++) {
    single << single_tcm[i].first << "\t" << (double)single_tcm[i].second/(double)single_tca[i].second << endl;
  }
  for (size_t i = 0; i < multi_tcm.size(); i++) {
    multi << multi_tcm[i].first << "\t" << (double)multi_tcm[i].second/(double)multi_tca[i].second << endl;
  }
  single.close();
  multi.close();
}

void test_BR_MSP() {

  int events[2] = {PAPI_BR_MSP, PAPI_BR_CN};
  ppvuill res = test_papi_events(events);

  ofstream single,multi;
  single.open("radix_single_BR_MSP.dat");
  multi.open("radix_multi_BR_MSP.dat");
  single << "#x\ty" << endl;
  multi << "#x\ty" << endl;
  vuill single_br_msp = res.first.first;
  vuill single_br_cn = res.first.second;
  vuill multi_br_msp = res.second.first;
  vuill multi_br_cn = res.second.second;
  //single
  for (size_t i = 0; i < single_br_msp.size(); i++) {
    single << single_br_msp[i].first << "\t" << (double)single_br_msp[i].second/(double)single_br_cn[i].second << endl;
  }
  for (size_t i = 0; i < multi_br_msp.size(); i++) {
    multi << multi_br_msp[i].first << "\t" << (double)multi_br_msp[i].second/(double)multi_br_cn[i].second << endl;
  }  single.close();
  multi.close();
}

void test_TLB() {

  int events[2] = {PAPI_TLB_DM, PAPI_TOT_INS};
  ppvuill res = test_papi_events(events);

  ofstream single,multi;
  single.open("radix_single_TLB_DM.dat");
  multi.open("radix_multi_TLB_DM.dat");
  single << "#x\ty" << endl;
  multi << "#x\ty" << endl;
  vuill single_TLB = res.first.first;
  vuill single_TOT = res.first.second;
  vuill multi_TLB = res.second.first;
  vuill multi_TOT = res.second.second;
  //single
  for (size_t i = 0; i < single_TLB.size(); i++) {
    single << single_TLB[i].first << "\t" << single_TLB[i].second << endl;
  }
  for (size_t i = 0; i < multi_TLB.size(); i++) {
    multi << multi_TLB[i].first << "\t" << multi_TLB[i].second << endl;
  }
  single.close();
  multi.close();
}

int main() {
  // srand(time(NULL));
  // test_TLB();
  // test_BR_MSP();
  // test_L2();
  // test_L3();
  // test_running_time();

  // return 0;

  size_t size = 100000000;
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
  radix_sort(data2, 4);  
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
