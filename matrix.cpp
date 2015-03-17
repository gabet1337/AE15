#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <stdio.h>
#include <algorithm>
#include <papi.h>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <time.h>
#include <string.h>

using namespace std;
using namespace std::chrono;

typedef vector<int> vi;
typedef vector< vi > matrix;
typedef pair<int,int> ii;

#define MAX_NUM 10

void resize_matrix(matrix &m, int rows, int cols) {
  m.assign(rows,vi(cols,0));
}

ii matrix_size(matrix &m) {
  return ii(m.size(),m[0].size());
}


void print_matrix(matrix (&m)) {
  size_t rows = m.size();
  size_t cols = m[0].size();
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      cout << m[i][j] << '\t';
    }
    cout << endl;
  }
}

void print_array(vector<int> a) {
  printf("size %d \n", (int)a.size());
  for (size_t i = 0; i < a.size(); i++)
    printf("%d ", a[i]);
}

void print_array_as_matrix(vector<int> a, int rows, int cols) {
  for (int i=0; i < rows; i++) {
    for (int j=0; j < cols; j++)
      cout << a[(j*rows)+i] << '\t';
    cout << endl;
  };
}

void preprocess_matrix_row_layout(matrix &mat, vector<int> &res) {
  size_t rows = matrix_size(mat).first;
  size_t cols = matrix_size(mat).second;
  res.resize(cols*rows,0);
  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < cols; j++)
      res[(i*cols)+j] = mat[i][j];
}

void preprocess_matrix_col_layout(matrix &mat, vector<int> &res) {
  size_t rows = matrix_size(mat).first;
  size_t cols = matrix_size(mat).second;
  res.resize(cols*rows,0);
  for (size_t i=0; i<cols;i++)
    for (size_t j=0; j<rows;j++)
      res[(i*rows)+j] = mat[j][i];
}

pair<vi,pair<vi,vi> > preprocess_mult_row(matrix &A, matrix &B) {
  
  vi a, b, c;

  int a_num_rows = matrix_size(A).first;
  int b_num_cols = matrix_size(B).second;

  c.resize(a_num_rows*b_num_cols,0);
  
  preprocess_matrix_row_layout(A,a);
  preprocess_matrix_col_layout(B,b);

  return make_pair(c, make_pair(a, b));
  
}

vector<int> mult_row(matrix &A, matrix &B, pair<vi,pair<vi,vi> > &prep_vec) {

  int a_num_rows = matrix_size(A).first;
  int a_num_cols = matrix_size(A).second;
  int b_num_rows = matrix_size(B).first;
  int b_num_cols = matrix_size(B).second;
  
  vi a = prep_vec.second.first;
  vi b = prep_vec.second.second;
  vi res = prep_vec.first;

  for (int k=0; k<b_num_cols;k++)
    for (int j=0; j<a_num_rows;j++)
      for (int i=0; i<a_num_cols;i++)
	res[(k*a_num_rows)+j] += (a[(j*a_num_cols)+i] * b[(k*b_num_rows)+i]);

  return res;
  
}

vector<int> mult_row_index_optimized(matrix &A, matrix &B, pair<vi,pair<vi,vi> > &prep_vec) {

  int a_num_rows = matrix_size(A).first;
  int a_num_cols = matrix_size(A).second;
  int b_num_cols = matrix_size(B).second;
  
  vi a = prep_vec.second.first;
  vi b = prep_vec.second.second;
  vi res = prep_vec.first;
  
  for (int k=0; k<b_num_cols;k++) {
    int m_offset = k*a_num_rows;
    for (int j=0; j<a_num_rows;j++) {
      int n_offset = j*a_num_cols;
      for (int i=0; i<a_num_cols;i++)
	res[m_offset+j] += (a[n_offset+i] * b[n_offset+i]);
    }
  }

  return res;
  
}

vector<int> mult_row_parallel(matrix &A, matrix &B, pair<vi,pair<vi,vi> > &prep_vec) {

  size_t a_num_rows = matrix_size(A).first;
  size_t a_num_cols = matrix_size(A).second;
  size_t b_num_cols = matrix_size(B).second;
  
  vi a = prep_vec.second.first;
  vi b = prep_vec.second.second;
  vi res = prep_vec.first;

  size_t i, j, k;
  
  //#pragma omp parallel for private(i,j,k)
  //#pragma omp parallel for ordered schedule(static)
#pragma omp parallel for default(none) shared(i,j,k, a_num_rows, b_num_cols, res, b, a, a_num_cols)
  for (k=0; k<b_num_cols;k++) {
    size_t m_offset = k*a_num_rows;
    for (j=0; j<a_num_rows;j++) {
      size_t n_offset = j*a_num_cols;
      for (i=0; i<a_num_cols;i++)
	res[m_offset+j] += (a[n_offset+i] * b[n_offset+i]);
    }
  }

  return res;
  
}

pair<vi,pair<vi,vi> > preprocess_mult_col(matrix &A, matrix &B) {
  
  vi a, b, c;

  int a_num_rows = matrix_size(A).first;
  int b_num_cols = matrix_size(B).second;

  c.resize(a_num_rows*b_num_cols,0);

  preprocess_matrix_col_layout(A,a);
  preprocess_matrix_row_layout(B,b);

  return make_pair(c, make_pair(a, b));
  
}

vector<int> mult_col(matrix &A, matrix &B, pair<vi,pair<vi,vi> > &prep_vec) {

  int a_num_rows = matrix_size(A).first;
  int a_num_cols = matrix_size(A).second;
  int b_num_cols = matrix_size(B).second;

  vi a = prep_vec.second.first;
  vi b = prep_vec.second.second;
  vi res = prep_vec.first;
  
  for (int k=0; k<b_num_cols;k++)
    for (int j=0; j<a_num_rows;j++)
      for (int i=0; i<a_num_cols;i++)
	res[(k*a_num_rows)+j] += (a[(i*a_num_rows)+j] * b[(i*b_num_cols)+k]);
 
  return res;
  
}

pair<vi,pair<vi,vi> > preprocess_mult_naive(matrix &A, matrix &B) {
  
  vi a, b, c;

  int a_num_rows = matrix_size(A).first;
  int b_num_cols = matrix_size(B).second;

  c.resize(a_num_rows*b_num_cols,0);
  
  preprocess_matrix_row_layout(A,a);
  preprocess_matrix_row_layout(B,b); 

  return make_pair(c, make_pair(a, b));
  
}





vector<int> mult_naive(matrix &A, matrix &B, pair<vi,pair<vi,vi> > &prep_vec) {

  int a_num_rows = matrix_size(A).first;
  int a_num_cols = matrix_size(A).second;
  int b_num_cols = matrix_size(B).second;

  vi a = prep_vec.second.first;
  vi b = prep_vec.second.second;
  vi res = prep_vec.first;
  
  for (int k=0; k<b_num_cols;k++)
    for (int j=0; j<a_num_rows;j++)
      for (int i=0; i<a_num_cols;i++)
	res[(k*a_num_rows)+j] = res[(k*a_num_rows)+j] + (a[(j*a_num_cols)+i] * b[(i*b_num_cols)+k]);
  
  return res;
  
}

vector<int> mult_naive_parallel(matrix &A, matrix &B, pair<vi,pair<vi,vi> > &prep_vec) {

  int a_num_rows = matrix_size(A).first;
  int a_num_cols = matrix_size(A).second;
  int b_num_cols = matrix_size(B).second;

  vi a = prep_vec.second.first;
  vi b = prep_vec.second.second;
  vi res = prep_vec.first;

  int k,j,i;
  //#pragma omp parallel for ordered schedule(static)
  //#pragma omp parallel for default(none) shared(i,j,k, a_num_rows, b_num_cols, res, b, a, a_num_cols)
#pragma omp parallel for private(i,j,k)
  for (k=0; k<b_num_cols;k++)
    for (j=0; j<a_num_rows;j++)
      for (i=0; i<a_num_cols;i++)
	res[(k*a_num_rows)+j] = res[(k*a_num_rows)+j] + (a[(j*a_num_cols)+i] * b[(i*b_num_cols)+k]);
  
  return res;
  
}

long long array_sum(vector<int> &a) {
  long long res = 0;
  for (size_t i=0; i<a.size(); i++)
    res += a[i];
  
  return res;
}

double test_runtime_mult_naive2(int size) {

  int** a = new int*[size];
  for(int i = 0; i < size; ++i)
    a[i] = new int[size];

  int** b = new int*[size];
  for(int i = 0; i < size; ++i)
    b[i] = new int[size];

  int** c = new int*[size];
  for(int i = 0; i < size; ++i)
    c[i] = new int[size];
  
  srand (time(NULL));
  
  // Initialize buffers.
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      int r = rand() % MAX_NUM;
      a[i][j] = r; //i*size+j+1; //  (float)i + j;
      r = rand() % MAX_NUM;
      b[i][j] = r;// i*size+j+1; //(float)i - j;
      c[i][j] = 0;
    }
  }

  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      for (int k = 0; k < size; ++k) {
	c[i][j] += a[i][k] * b[k][j];
      }
    }
  }

  high_resolution_clock::time_point t2 = high_resolution_clock::now();

  return (double) std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
  
}

double test_runtime_mult_parallel_naive2(int size) {

  int** a = new int*[size];
  for(int i = 0; i < size; ++i)
    a[i] = new int[size];

  int** b = new int*[size];
  for(int i = 0; i < size; ++i)
    b[i] = new int[size];

  int** c = new int*[size];
  for(int i = 0; i < size; ++i)
    c[i] = new int[size];
  
  srand (time(NULL));
  
  // Initialize buffers.
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      int r = rand() % MAX_NUM;
      a[i][j] = r; //i*size+j+1; //  (float)i + j;
      r = rand() % MAX_NUM;
      b[i][j] = r;// i*size+j+1; //(float)i - j;
      c[i][j] = 0;
    }
  }


  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
#pragma omp parallel for default(none) shared(a,b,c,size)
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      for (int k = 0; k < size; ++k) {
	c[i][j] += a[i][k] * b[k][j];
      }
    }
  }

  high_resolution_clock::time_point t2 = high_resolution_clock::now();

  return (double) std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
  
}

matrix sum_matrix(matrix A, matrix B) {
  
  matrix res;
  resize_matrix(res, matrix_size(A).first, matrix_size(A).second);
  
  for (int i=0; i < matrix_size(A).first; i++)
    for (int j=0; j < matrix_size(A).second; j++)
      res[i][j] = A[i][j] + B[i][j];
  
  return res;

}

matrix hor_concat_matrix(matrix &A, matrix &B) {

  matrix res;
  resize_matrix(res, matrix_size(A).first+matrix_size(B).first, matrix_size(A).second);
  
  for (int i=0; i < matrix_size(A).first; i++){
    res[i] = A[i];
  }
  
  for (int i=0; i < matrix_size(B).first; i++){
    res[i+matrix_size(A).first] = B[i];
  }
  
  return res;
  
}

matrix ver_concat_matrix(matrix &A, matrix &B) {

  matrix res;
  resize_matrix(res, matrix_size(A).first, matrix_size(A).second+matrix_size(B).second);
  
  for (int i=0; i<matrix_size(A).first; i++) {
    for (int j=0; j<matrix_size(A).second; j++) {
      res[i][j] = A[i][j];
    }
  }
  
  for (int i=0; i<matrix_size(B).first; i++) {
    for (int j=matrix_size(B).second; j<matrix_size(A).second+matrix_size(B).second; j++) {
      res[i][j] = B[i][j-matrix_size(B).second];
    }
  }
  
  return res;
}

matrix hor_split_matrix(matrix &M, int start, int end) {
  
  matrix res;
  resize_matrix(res, end-start+1, matrix_size(M).second);
  
  for (int i=start; i < end+1; i++){
    res[i-start] = M[i];
  }
  
  return res;
  
}

matrix ver_split_matrix(matrix &M, int start, int end) {
  
  matrix res;
  resize_matrix(res, matrix_size(M).first, end-start+1);
  
  for (int i=0; i<matrix_size(M).first; i++) {
    for (int j=start; j<end+1; j++) {
      res[i][j-start] = M[i][j];
    }
  }
  
  return res;
  
}

void print_to_plot(vector<pair<int,double> > &data, const char* fname) {
  ofstream myfile;
  myfile.open(fname);
  //cout << "#x\ty" << endl;
  myfile << "#x\ty" << endl;
  for (size_t i = 0; i < data.size(); i++) {
    myfile << data[i].first << "\t" << data[i].second << endl;
    //cout << data[i].first << "\t" << data[i].second << endl;
  }
  myfile.close();
}

matrix mult_rec(matrix &A, matrix &B, matrix C, int m_start, int m_end, int n_start, int n_end, int p_start,  int p_end) {
  
  // A_mxn, B_nxp, C_mxp
  
  int m = m_end - m_start + 1;
  int n = n_end - n_start + 1;
  int p = p_end - p_start + 1;
  
  if ((m == 1) && (n == 1) && (p == 1)) {
    C[0][0] = C[0][0] + (A[m_start-1][n_start-1] * B[n_start-1][p_start-1]);
    return C;
  } else if (m >= max(n, p)) {
    
    int m_mid = m/2; // split m
    
    matrix C1 = hor_split_matrix(C, 0, m_mid-1); // Copy upper half of C to C1
    matrix a1b = mult_rec(A, B, C1, m_start, m_end-m_mid, n_start, n_end, p_start, p_end); // result of this is A_1 * B
    
    matrix C2 = hor_split_matrix(C, m_mid, m-1);
    matrix a2b = mult_rec(A, B, C2, m_end-m_mid+1, m_end, n_start, n_end, p_start, p_end);
    
    C1 = sum_matrix(C1, a1b);
    C2 = sum_matrix(C2, a2b);
    
    return hor_concat_matrix(C1, C2);
    
  } else if (n >= max(m, p)) {
    
    int n_mid = n/2; // split n
    
    return sum_matrix(mult_rec(A, B, C, m_start, m_end, n_start, n_end-n_mid, p_start, p_end), mult_rec(A, B, C, m_start, m_end, n_end-n_mid+1, n_end, p_start, p_end));
    
  } else if (p >= max(m, n)) {
    
    int p_mid = p/2; //split p
    
    matrix C1 = ver_split_matrix(C, 0, p_mid-1);
    matrix ab1 = mult_rec(A, B, C1, m_start, m_end, n_start, n_end, p_start, p_end-p_mid);
    
    matrix C2 = ver_split_matrix(C, p_mid, p-1);
    matrix ab2 = mult_rec(A, B, C1, m_start, m_end, n_start, n_end, p_end-p_mid+1, p_end);
    
    C1 = sum_matrix(C1, ab1);
    C2 = sum_matrix(C2, ab2);
    
    return ver_concat_matrix(C1, C2);

  }
  return C;
}

matrix mult_rec_cutoff(matrix &A, matrix &B, matrix C, int m_start, int m_end, int n_start, int n_end, int p_start,  int p_end, int cutoff) {
  
  // A_mxn, B_nxp, C_mxp
  
  int m = m_end - m_start + 1;
  int n = n_end - n_start + 1;
  int p = p_end - p_start + 1;
  // cout << m_start << " " << m_end << " " << n_start << " " << n_end << " " << p_start << " " << p_end << endl;
  if ((m <= cutoff) && (n <= cutoff) && (p <= cutoff)) {
    m_start--; n_start--; p_start--;
    // cout << m_start << " " << m_end << " " << n_start << " " << n_end << " " << p_start << " " << p_end << endl;
    for (int i = m_start; i < m_end; ++i)
      for (int j = n_start; j < n_end; ++j)
	for (int k = p_start; k < p_end; ++k)
	  C[i-m_start][j-n_start] += A[i][k] * B[k][j];
  } else if (m >= max(n, p)) {
    
    int m_mid = m/2; // split m
    
    matrix C1 = hor_split_matrix(C, 0, m_mid-1); // Copy upper half of C to C1
    matrix a1b = mult_rec_cutoff(A, B, C1, m_start, m_end-m_mid, n_start, n_end, p_start, p_end,cutoff); // result of this is A_1 * B
    
    matrix C2 = hor_split_matrix(C, m_mid, m-1);
    matrix a2b = mult_rec_cutoff(A, B, C2, m_end-m_mid+1, m_end, n_start, n_end, p_start, p_end,cutoff);
    
    C1 = sum_matrix(C1, a1b);
    C2 = sum_matrix(C2, a2b);
    
    return hor_concat_matrix(C1, C2);
    
  } else if (n >= max(m, p)) {
    
    int n_mid = n/2; // split n
    
    return sum_matrix(mult_rec_cutoff(A, B, C, m_start, m_end, n_start, n_end-n_mid, p_start, p_end,cutoff), mult_rec_cutoff(A, B, C, m_start, m_end, n_end-n_mid+1, n_end, p_start, p_end,cutoff));
    
  } else if (p >= max(m, n)) {
    
    int p_mid = p/2; //split p
    
    matrix C1 = ver_split_matrix(C, 0, p_mid-1);
    matrix ab1 = mult_rec_cutoff(A, B, C1, m_start, m_end, n_start, n_end, p_start, p_end-p_mid,cutoff);
    
    matrix C2 = ver_split_matrix(C, p_mid, p-1);
    matrix ab2 = mult_rec_cutoff(A, B, C1, m_start, m_end, n_start, n_end, p_end-p_mid+1, p_end,cutoff);
    
    C1 = sum_matrix(C1, ab1);
    C2 = sum_matrix(C2, ab2);
    
    return ver_concat_matrix(C1, C2);

  }
  return C;
}

pair<double, pair<double,double> > test_mult_rec(matrix &A, matrix &B, int events[], int event_size, int numruns) {

  matrix C;
  resize_matrix(C, matrix_size(A).first, matrix_size(B).second);
  
  // Clock, PAPI
  clock_t start,end;
  double timer = 0;
  
  PAPI_library_init(PAPI_VER_CURRENT);
  long long values[event_size];

  long long count0 = 0;
  long long count1 = 0;
  
  for (int i=0; i<numruns; i++) {

    srand (time(NULL));

    for (int i=0; i<matrix_size(A).first; i++)
      for (int j=0; j<matrix_size(A).second; j++) {
	int r = rand() % MAX_NUM;
	A[i][j] = r; // (i*matrix_size(A).second)+j+1;
      }

    for (int i=0; i<matrix_size(B).first; i++)
      for (int j=0; j<matrix_size(B).second; j++) {
	int r = rand() % MAX_NUM;
	B[i][j] = r; // (i*matrix_size(B).second)+j+1;
      }
    
    PAPI_start_counters(events, event_size);

    start = clock();

    // Multiplying matrices A, B
    matrix res = mult_rec(A, B, C, 1, matrix_size(A).first, 1, matrix_size(B).first, 1, matrix_size(B).second);
  
    end = clock();

    PAPI_stop_counters(values, event_size);
  
    count0 += values[0];
  
    if (event_size == 2)
      count1 += values[1];
  
    timer += (end - start) / (double)(CLOCKS_PER_SEC / 1000);

  }
  
  return make_pair(timer/numruns, make_pair(count0/numruns, count1/numruns));
  
  //print_matrix(res);

}


pair<double, pair<double,double> > test_mult_rec_cutoff(matrix &A, matrix &B, int events[], int event_size, int numruns, int cutoff) {

  matrix C;
  resize_matrix(C, matrix_size(A).first, matrix_size(B).second);
  
  // Clock, PAPI
  clock_t start,end;
  double timer = 0;
  
  PAPI_library_init(PAPI_VER_CURRENT);
  long long values[event_size];

  long long count0 = 0;
  long long count1 = 0;
  
  for (int i=0; i<numruns; i++) {

    srand (time(NULL));

    for (int i=0; i<matrix_size(A).first; i++)
      for (int j=0; j<matrix_size(A).second; j++) {
	int r = rand() % MAX_NUM;
	A[i][j] = r; // (i*matrix_size(A).second)+j+1;
      }

    for (int i=0; i<matrix_size(B).first; i++)
      for (int j=0; j<matrix_size(B).second; j++) {
	int r = rand() % MAX_NUM;
	B[i][j] = r; // (i*matrix_size(B).second)+j+1;
      }
    
    PAPI_start_counters(events, event_size);

    start = clock();

    // Multiplying matrices A, B
    matrix res = mult_rec_cutoff(A, B, C, 1, matrix_size(A).first, 1, matrix_size(B).first, 1, matrix_size(B).second, cutoff);
  
    end = clock();

    PAPI_stop_counters(values, event_size);
  
    count0 += values[0];
  
    if (event_size == 2)
      count1 += values[1];
  
    timer += (end - start) / (double)(CLOCKS_PER_SEC / 1000);

  }
  
  return make_pair(timer/numruns, make_pair(count0/numruns, count1/numruns));
  
  //print_matrix(res);

}




pair<double, pair<double,double> > test_row_mult_exclude_preprocess(matrix &A, matrix &B, int events[], int event_size, int numruns) {
  
  // Preprocess matrices A, B
  pair<vi,pair<vi,vi> > prv = preprocess_mult_row(A, B);

  // Clock, PAPI
  double timer = 0;
  
  PAPI_library_init(PAPI_VER_CURRENT);
  long long values[event_size];

  long long count0 = 0;
  long long count1 = 0;
  
  for (int i=0; i<numruns; i++) {

    srand (time(NULL));

    for (int i=0; i<matrix_size(A).first; i++)
      for (int j=0; j<matrix_size(A).second; j++) {
	int r = rand() % MAX_NUM;
	A[i][j] = r; // (i*matrix_size(A).second)+j+1;
      }

    for (int i=0; i<matrix_size(B).first; i++)
      for (int j=0; j<matrix_size(B).second; j++) {
	int r = rand() % MAX_NUM;
	B[i][j] = r; // (i*matrix_size(B).second)+j+1;
      }
    
    PAPI_start_counters(events, event_size);

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    // Multiplying matrices A, B
    vector<int> res_row = mult_row(A, B, prv);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    PAPI_stop_counters(values, event_size);
  
    count0 += values[0];
  
    if (event_size == 2)
      count1 += values[1];

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
    
    timer += (double) duration / 1000;

  }
  
  return make_pair(timer/numruns, make_pair(count0/numruns, count1/numruns));
    
}

pair<double, pair<double,double> > test_row_mult_idx_opt_exclude_preprocess(matrix &A, matrix &B, int events[], int event_size, int numruns) {
  
  // Preprocess matrices A, B
  pair<vi,pair<vi,vi> > prv = preprocess_mult_row(A, B);

  // Clock, PAPI
  clock_t start,end;
  double timer = 0;
  
  PAPI_library_init(PAPI_VER_CURRENT);
  long long values[event_size];

  long long count0 = 0;
  long long count1 = 0;
  
  for (int i=0; i<numruns; i++) {

    srand (time(NULL));

    for (int i=0; i<matrix_size(A).first; i++)
      for (int j=0; j<matrix_size(A).second; j++) {
	int r = rand() % MAX_NUM;
	A[i][j] = r; // (i*matrix_size(A).second)+j+1;
      }

    for (int i=0; i<matrix_size(B).first; i++)
      for (int j=0; j<matrix_size(B).second; j++) {
	int r = rand() % MAX_NUM;
	B[i][j] = r; // (i*matrix_size(B).second)+j+1;
      }
    
    PAPI_start_counters(events, event_size);

    start = clock();

    // Multiplying matrices A, B
    vector<int> res_row = mult_row_index_optimized(A, B, prv);
  
    end = clock();

    PAPI_stop_counters(values, event_size);
  
    count0 += values[0];
  
    if (event_size == 2)
      count1 += values[1];
  
    timer += (end - start) / (double)(CLOCKS_PER_SEC / 1000);

  }
  
  return make_pair(timer/numruns, make_pair(count0/numruns, count1/numruns));
    
}

pair<double, pair<double,double> > test_row_mult_parallel_exclude_preprocess(matrix &A, matrix &B, int events[], int event_size, int numruns) {
  
  // Preprocess matrices A, B
  pair<vi,pair<vi,vi> > prv = preprocess_mult_row(A, B);

  // Clock, PAPI
  clock_t start,end;
  double timer = 0;
  
  PAPI_library_init(PAPI_VER_CURRENT);
  long long values[event_size];

  long long count0 = 0;
  long long count1 = 0;
  
  for (int i=0; i<numruns; i++) {

    srand (time(NULL));

    for (int i=0; i<matrix_size(A).first; i++)
      for (int j=0; j<matrix_size(A).second; j++) {
	int r = rand() % MAX_NUM;
	A[i][j] = r; // (i*matrix_size(A).second)+j+1;
      }

    for (int i=0; i<matrix_size(B).first; i++)
      for (int j=0; j<matrix_size(B).second; j++) {
	int r = rand() % MAX_NUM;
	B[i][j] = r; // (i*matrix_size(B).second)+j+1;
      }
    
    PAPI_start_counters(events, event_size);

    start = clock();

    // Multiplying matrices A, B
    vector<int> res_row = mult_row_parallel(A, B, prv);
  
    end = clock();

    PAPI_stop_counters(values, event_size);
  
    count0 += values[0];
  
    if (event_size == 2)
      count1 += values[1];
  
    timer += (end - start) / (double)(CLOCKS_PER_SEC / 1000);

  }
  
  return make_pair(timer/numruns, make_pair(count0/numruns, count1/numruns));
    
}

pair<double, pair<double,double> > test_col_mult_exclude_preprocess(matrix &A, matrix &B, int events[], int event_size, int numruns) {
  
  // Preprocess matrices A, B
  pair<vi,pair<vi,vi> > prv = preprocess_mult_col(A, B);

  // Clock, PAPI
  double timer = 0;
  
  PAPI_library_init(PAPI_VER_CURRENT);
  long long values[event_size];

  long long count0 = 0;
  long long count1 = 0;
  
  for (int i=0; i<numruns; i++) {

    srand (time(NULL));

    for (int i=0; i<matrix_size(A).first; i++)
      for (int j=0; j<matrix_size(A).second; j++) {
	int r = rand() % MAX_NUM;
	A[i][j] = r; // (i*matrix_size(A).second)+j+1;
      }

    for (int i=0; i<matrix_size(B).first; i++)
      for (int j=0; j<matrix_size(B).second; j++) {
	int r = rand() % MAX_NUM;
	B[i][j] = r; // (i*matrix_size(B).second)+j+1;
      }
    
    PAPI_start_counters(events, event_size);

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    // Multiplying matrices A, B
    vector<int> res = mult_col(A, B, prv);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    PAPI_stop_counters(values, event_size);

    count0 += values[0];
  
    if (event_size == 2)
      count1 += values[1];

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
    
    timer += (double) duration / 1000;

    //print_array_as_matrix(res, matrix_size(A).first, matrix_size(B).second);

  }
  
  return make_pair(timer/numruns, make_pair(count0/numruns, count1/numruns));
    
}

/*
      for (int i=1; i<10; i++) {
      int size = i*200;
      cout << "N = " << size << " ----------------------------" << endl;
      cout << "Single:\t\t" << test_runtime_mult_naive2(i*200) / 1000 << endl;
      omp_set_num_threads(2);
      cout << "Par (2 cores):\t" << test_runtime_mult_parallel_naive2(i*200) / 1000 << endl;
      omp_set_num_threads(4);
      cout << "Par (4 cores):\t" << test_runtime_mult_parallel_naive2(i*200) / 1000 << endl;
      cout << endl;
      }
     */

pair<double, pair<double,double> > test_naive2_mult_exclude_preprocess(int size, int events[], int event_size, int numruns) {
  
  // Clock, PAPI
  double timer = 0;
  
  PAPI_library_init(PAPI_VER_CURRENT);
  long long values[event_size];

  long long count0 = 0;
  long long count1 = 0;
  
  for (int i=0; i<numruns; i++) {

    int** a = new int*[size];
    for(int i = 0; i < size; ++i)
      a[i] = new int[size];
    
    int** b = new int*[size];
    for(int i = 0; i < size; ++i)
      b[i] = new int[size];
    
    int** c = new int*[size];
    for(int i = 0; i < size; ++i)
      c[i] = new int[size];
  
    srand (time(NULL));
  
    // Initialize buffers.
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
	int r = rand() % MAX_NUM;
	a[i][j] = r; //i*size+j+1; //  (float)i + j;
	r = rand() % MAX_NUM;
	b[i][j] = r;// i*size+j+1; //(float)i - j;
	c[i][j] = 0;
      }
    }
    
    PAPI_start_counters(events, event_size);
    
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    
    // Multiplying matrices A, B
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
	for (int k = 0; k < size; ++k) {
	  c[i][j] += a[i][k] * b[k][j];
	}
      }
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    PAPI_stop_counters(values, event_size);

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
    
    count0 += values[0];
  
    if (event_size == 2)
      count1 += values[1];
  
    timer += (double) duration / 1000;

    //print_array_as_matrix(res, matrix_size(A).first, matrix_size(B).second);

  }
  
  return make_pair(timer/numruns, make_pair(count0/numruns, count1/numruns));
    
}

pair<double, pair<double,double> > test_naive2_parallel_mult_exclude_preprocess(int size, int numthreads, int events[], int event_size, int numruns) {

  omp_set_num_threads(numthreads);
  
  // Clock, PAPI
  double timer = 0;
  
  PAPI_library_init(PAPI_VER_CURRENT);
  long long values[event_size];

  long long count0 = 0;
  long long count1 = 0;
  
  for (int i=0; i<numruns; i++) {

    int** a = new int*[size];
    for(int i = 0; i < size; ++i)
      a[i] = new int[size];
    
    int** b = new int*[size];
    for(int i = 0; i < size; ++i)
      b[i] = new int[size];
    
    int** c = new int*[size];
    for(int i = 0; i < size; ++i)
      c[i] = new int[size];
  
    srand (time(NULL));
  
    // Initialize buffers.
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
	int r = rand() % MAX_NUM;
	a[i][j] = r; //i*size+j+1; //  (float)i + j;
	r = rand() % MAX_NUM;
	b[i][j] = r;// i*size+j+1; //(float)i - j;
	c[i][j] = 0;
      }
    }
    
    PAPI_start_counters(events, event_size);
    
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    
    // Multiplying matrices A, B
#pragma omp parallel for default(none) shared(a,b,c,size)
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
	for (int k = 0; k < size; ++k) {
	  c[i][j] += a[i][k] * b[k][j];
	}
      }
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    PAPI_stop_counters(values, event_size);

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
    
    count0 += values[0];
  
    if (event_size == 2)
      count1 += values[1];
  
    timer += (double) duration / 1000;

    //print_array_as_matrix(res, matrix_size(A).first, matrix_size(B).second);

  }
  
  return make_pair(timer/numruns, make_pair(count0/numruns, count1/numruns));
    
}

pair<double, pair<double,double> > test_naive_mult_exclude_preprocess(matrix &A, matrix &B, int events[], int event_size, int numruns) {
  
  // Preprocess matrices A, B
  pair<vi,pair<vi,vi> > prv = preprocess_mult_naive(A, B);

  // Clock, PAPI
  double timer = 0;
  
  PAPI_library_init(PAPI_VER_CURRENT);
  long long values[event_size];

  long long count0 = 0;
  long long count1 = 0;
  
  for (int i=0; i<numruns; i++) {

    srand (time(NULL));

    for (int i=0; i<matrix_size(A).first; i++)
      for (int j=0; j<matrix_size(A).second; j++) {
	int r = rand() % MAX_NUM;
	A[i][j] = r; // (i*matrix_size(A).second)+j+1;
      }

    for (int i=0; i<matrix_size(B).first; i++)
      for (int j=0; j<matrix_size(B).second; j++) {
	int r = rand() % MAX_NUM;
	B[i][j] = r; // (i*matrix_size(B).second)+j+1;
      }
    
    PAPI_start_counters(events, event_size);
    
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    
    // Multiplying matrices A, B
    vector<int> res = mult_naive(A, B, prv);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    PAPI_stop_counters(values, event_size);

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
    
    count0 += values[0];
  
    if (event_size == 2)
      count1 += values[1];
  
    timer += (double) duration / 1000;

    //print_array_as_matrix(res, matrix_size(A).first, matrix_size(B).second);

  }
  
  return make_pair(timer/numruns, make_pair(count0/numruns, count1/numruns));
    
}

pair<double, pair<double,double> > test_naive_mult_parallel_exclude_preprocess(matrix &A, matrix &B, int events[], int event_size, int numruns) {
  
  // Preprocess matrices A, B
  pair<vi,pair<vi,vi> > prv = preprocess_mult_naive(A, B);

  // Clock, PAPI
  clock_t start,end;
  double timer = 0;
  
  PAPI_library_init(PAPI_VER_CURRENT);
  long long values[event_size];

  long long count0 = 0;
  long long count1 = 0;
  
  for (int i=0; i<numruns; i++) {

    srand (time(NULL));

    for (int i=0; i<matrix_size(A).first; i++)
      for (int j=0; j<matrix_size(A).second; j++) {
	int r = rand() % MAX_NUM;
	A[i][j] = r; // (i*matrix_size(A).second)+j+1;
      }

    for (int i=0; i<matrix_size(B).first; i++)
      for (int j=0; j<matrix_size(B).second; j++) {
	int r = rand() % MAX_NUM;
	B[i][j] = r; // (i*matrix_size(B).second)+j+1;
      }
    
    PAPI_start_counters(events, event_size);

    start = clock();

    // Multiplying matrices A, B
    vector<int> res = mult_naive_parallel(A, B, prv);
  
    end = clock();

    PAPI_stop_counters(values, event_size);
  
    count0 += values[0];
  
    if (event_size == 2)
      count1 += values[1];
  
    timer += (end - start) / (double)(CLOCKS_PER_SEC / 1000);

    //print_array_as_matrix(res, matrix_size(A).first, matrix_size(B).second);

  }
  
  return make_pair(timer/numruns, make_pair(count0/numruns, count1/numruns));
    
}

void test_runtime_l2_l3_accesses() {

  vector<pair<int, double > > naive_runtime, naive2_runtime, naive2_par2_runtime, naive2_par4_runtime, col_runtime, row_runtime;
  vector<pair<int, double > > naive_l2, naive2_par2_l2, naive2_par4_l2, naive2_l2, col_l2, row_l2;
  vector<pair<int, double > > naive_l3, naive2_par2_l3, naive2_par4_l3, naive2_l3, col_l3, row_l3;

  int size;

  int testruns = 1;
  int events[2] = {PAPI_L2_TCA, PAPI_L3_TCA};
  int event_size = 2;

  for (int run=128; run<4097; run+=128) {

    size = run;

    cout << "Running test n = " << run << endl;
    
    int m = size;
    int n = size;
    int p = size;
  
    matrix A, B;
  
    resize_matrix(A, m, n);
    resize_matrix(B, n, p);

    cout << " - starting test of row_mult" << endl;
    pair<double,pair<double,double> > res_row = test_row_mult_exclude_preprocess(A, B, events, event_size, testruns);
    row_runtime.push_back(make_pair(run, res_row.first));
    row_l2.push_back(make_pair(run, res_row.second.first));
    row_l3.push_back(make_pair(run, res_row.second.second));

    if (run <= 1024) {
      cout << " - starting test of col_mult" << endl;
      pair<double,pair<double,double> > res_col = test_col_mult_exclude_preprocess(A, B, events, event_size, testruns);
      col_runtime.push_back(make_pair(run, res_col.first));
      col_l2.push_back(make_pair(run, res_col.second.first));
      col_l3.push_back(make_pair(run, res_col.second.second));
    }

    if (run <= 2048) {
    cout << " - starting test of naive_mult" << endl;
    pair<double,pair<double,double> > res_naive = test_naive_mult_exclude_preprocess(A, B, events, event_size, testruns);
    naive_runtime.push_back(make_pair(run, res_naive.first));
    naive_l2.push_back(make_pair(run, res_naive.second.first));
    naive_l3.push_back(make_pair(run, res_naive.second.second));    

    cout << " - starting test of naive2_mult" << endl;
    pair<double,pair<double,double> > res_naive2 = test_naive2_mult_exclude_preprocess(size, events, event_size, testruns);
    naive2_runtime.push_back(make_pair(run, res_naive2.first));
    naive2_l2.push_back(make_pair(run, res_naive2.second.first));
    naive2_l3.push_back(make_pair(run, res_naive2.second.second));    

    cout << " - starting test of naive2_par2_mult" << endl;
    pair<double,pair<double,double> > res_naive2_par2 = test_naive2_parallel_mult_exclude_preprocess(size, 2, events, event_size, testruns);
    naive2_par2_runtime.push_back(make_pair(run, res_naive2_par2.first));
    naive2_par2_l2.push_back(make_pair(run, res_naive2_par2.second.first));
    naive2_par2_l3.push_back(make_pair(run, res_naive2_par2.second.second));    

    cout << " - starting test of naive2_par4_mult" << endl;
    pair<double,pair<double,double> > res_naive2_par4 = test_naive2_parallel_mult_exclude_preprocess(size, 4, events, event_size, testruns);
    naive2_par4_runtime.push_back(make_pair(run, res_naive2_par4.first));
    naive2_par4_l2.push_back(make_pair(run, res_naive2_par4.second.first));
    naive2_par4_l3.push_back(make_pair(run, res_naive2_par4.second.second));
    }

  }

  print_to_plot(naive_runtime, "data_matrix/naive_runtime.dat");
  print_to_plot(naive2_runtime, "data_matrix/naive2_runtime.dat");
  print_to_plot(naive2_par2_runtime, "data_matrix/naive2_par2_runtime.dat");
  print_to_plot(naive2_par4_runtime, "data_matrix/naive2_par4_runtime.dat");
  print_to_plot(col_runtime, "data_matrix/col_runtime.dat");
  print_to_plot(row_runtime, "data_matrix/row_runtime.dat");

  print_to_plot(naive_l2, "data_matrix/naive_l2.dat");
  print_to_plot(naive2_l2, "data_matrix/naive2_l2.dat");
  print_to_plot(naive2_par2_l2, "data_matrix/naive2_par2_l2.dat");
  print_to_plot(naive2_par4_l2, "data_matrix/naive_par4_l2.dat");
  print_to_plot(col_l2, "data_matrix/col_l2.dat");
  print_to_plot(row_l2, "data_matrix/row_l2.dat");

  print_to_plot(naive_l3, "data_matrix/naive_l3.dat");
  print_to_plot(naive2_l3, "data_matrix/naive2_l3.dat");
  print_to_plot(naive2_par2_l3, "data_matrix/naive_par2_l3.dat");
  print_to_plot(naive2_par4_l3, "data_matrix/naive_par4_l3.dat");
  print_to_plot(col_l3, "data_matrix/col_l3.dat");
  print_to_plot(row_l3, "data_matrix/row_l3.dat");

}

void test_runtime_mirror(int from, int to, int testruns, int events[], int event_size, const char* fname, const char* title, const char* count0_title, const char* count1_title) {

  vector<pair<int, double > > naive_runtime, naive2_runtime, naive2_par2_runtime, naive2_par4_runtime, col_runtime, row_runtime;
  vector<pair<int, double > > naive_l2, naive2_par2_l2, naive2_par4_l2, naive2_l2, col_l2, row_l2;
  vector<pair<int, double > > naive_l3, naive2_par2_l3, naive2_par4_l3, naive2_l3, col_l3, row_l3;

  int size;

  //int testruns = 1;
  //int events[2] = {PAPI_L2_TCA, PAPI_L3_TCA};
  //int event_size = 2;

  for (int run=from; run<to+1; run++) {

    size = run;

    cout << "Running test n = " << run << endl;
    
    int m = run;
    int n = 4096/4;
    int p = run;
  
    matrix A, B;
  
    resize_matrix(A, m, n);
    resize_matrix(B, n, p);

    cout << " - starting test of row_mult" << endl;
    pair<double,pair<double,double> > res_row = test_row_mult_exclude_preprocess(A, B, events, event_size, testruns);
    row_runtime.push_back(make_pair(run, res_row.first));
    row_l2.push_back(make_pair(run, res_row.second.first));
    row_l3.push_back(make_pair(run, res_row.second.second));

    if (run <= 200) {
      cout << " - starting test of col_mult" << endl;
      pair<double,pair<double,double> > res_col = test_col_mult_exclude_preprocess(A, B, events, event_size, testruns);
      col_runtime.push_back(make_pair(run, res_col.first));
      col_l2.push_back(make_pair(run, res_col.second.first));
      col_l3.push_back(make_pair(run, res_col.second.second));
    }

    if (run <= 200) {
    cout << " - starting test of naive_mult" << endl;
    pair<double,pair<double,double> > res_naive = test_naive_mult_exclude_preprocess(A, B, events, event_size, testruns);
    naive_runtime.push_back(make_pair(run, res_naive.first));
    naive_l2.push_back(make_pair(run, res_naive.second.first));
    naive_l3.push_back(make_pair(run, res_naive.second.second));    
    }
    
    cout << " - starting test of naive2_mult" << endl;
    pair<double,pair<double,double> > res_naive2 = test_naive2_mult_exclude_preprocess(size, events, event_size, testruns);
    naive2_runtime.push_back(make_pair(run, res_naive2.first));
    naive2_l2.push_back(make_pair(run, res_naive2.second.first));
    naive2_l3.push_back(make_pair(run, res_naive2.second.second));    

    cout << " - starting test of naive2_par2_mult" << endl;
    pair<double,pair<double,double> > res_naive2_par2 = test_naive2_parallel_mult_exclude_preprocess(size, 2, events, event_size, testruns);
    naive2_par2_runtime.push_back(make_pair(run, res_naive2_par2.first));
    naive2_par2_l2.push_back(make_pair(run, res_naive2_par2.second.first));
    naive2_par2_l3.push_back(make_pair(run, res_naive2_par2.second.second));    

    cout << " - starting test of naive2_par4_mult" << endl;
    pair<double,pair<double,double> > res_naive2_par4 = test_naive2_parallel_mult_exclude_preprocess(size, 4, events, event_size, testruns);
    naive2_par4_runtime.push_back(make_pair(run, res_naive2_par4.first));
    naive2_par4_l2.push_back(make_pair(run, res_naive2_par4.second.first));
    naive2_par4_l3.push_back(make_pair(run, res_naive2_par4.second.second));

  }

  char fname1[100];
  strcpy(fname1,fname);
  strcat(fname1, "naive_");
  strcat(fname1, title);
  strcat(fname1, "_runtime.dat");
  print_to_plot(naive_runtime, fname1); // strcat(fname, "/1_naive_l2_l3_tca_runtime.dat"));

  char fname2[100];
  strcpy(fname2,fname);
  strcat(fname2, "naive2_");
  strcat(fname2, title);
  strcat(fname2, "_runtime.dat");
  print_to_plot(naive2_runtime, fname2); // fname + "/1_naive2_l2_l3_tca_runtime.dat");

  char fname3[100];
  strcpy(fname3,fname);
  strcat(fname3, "naive2_par2_");
  strcat(fname3, title);
  strcat(fname3, "_runtime.dat");
  print_to_plot(naive2_par2_runtime, fname3); // fname + "/1_naive2_par2_l2_l3_tca_runtime.dat");

  char fname4[100];
  strcpy(fname4,fname);
  strcat(fname4, "naive2_par4_");
  strcat(fname4, title);
  strcat(fname4, "_runtime.dat");
  print_to_plot(naive2_par4_runtime, fname4); // fname + "/1_naive2_par4_l2_l3_tca_runtime.dat");

  char fname5[100];
  strcpy(fname5,fname);
  strcat(fname5, "a-trans_");
  strcat(fname5, title);
  strcat(fname5, "_runtime.dat");
  print_to_plot(col_runtime, fname5); // fname + "/1_a-trans_l2_l3_tca_runtime.dat");

  char fname6[100];
  strcpy(fname6,fname);
  strcat(fname6, "b-trans_");
  strcat(fname6, title);
  strcat(fname6, "_runtime.dat");
  print_to_plot(row_runtime, fname6); // fname + "/1_b-trans_l2_l3_tca_runtime.dat");

  // PAPI Counter1
  char fname7[100];
  strcpy(fname7, fname);
  strcat(fname7, "naive_");
  strcat(fname7, title);
  strcat(fname7, "_");
  strcat(fname7, count0_title);
  strcat(fname7, ".dat");
  print_to_plot(naive_l2, fname7); // fname + "/2_naive_l2_tca.dat");

  char fname8[100];
  strcpy(fname8, fname);
  strcat(fname8, "naive2_");
  strcat(fname8, title);
  strcat(fname8, "_");
  strcat(fname8, count0_title);
  strcat(fname8, ".dat");
  print_to_plot(naive2_l2, fname8); // fname + "/2_naive2_l2_tca.dat");

  char fname9[100];
  strcpy(fname9, fname);
  strcat(fname9, "naive2_par2_");
  strcat(fname9, title);
  strcat(fname9, "_");
  strcat(fname9, count0_title);
  strcat(fname9, ".dat");
  print_to_plot(naive2_par2_l2, fname9); // fname + "/2_naive2_par2_l2_tca.dat");

  char fname10[100];
  strcpy(fname10, fname);
  strcat(fname10, "naive2_par4_");
  strcat(fname10, title);
  strcat(fname10, "_");
  strcat(fname10, count0_title);
  strcat(fname10, ".dat");
  print_to_plot(naive2_par4_l2, fname10); // fname + "/2_naive_par4_l2_tca.dat");

  char fname11[100];
  strcpy(fname11, fname);
  strcat(fname11, "a-trans_");
  strcat(fname11, title);
  strcat(fname11, "_");
  strcat(fname11, count0_title);
  strcat(fname11, ".dat");
  print_to_plot(col_l2, fname11); // fname + "/2_a-trans_l2_tca.dat");

  char fname12[100];
  strcpy(fname12, fname);
  strcat(fname12, "b-trans_");
  strcat(fname12, title);
  strcat(fname12, "_");
  strcat(fname12, count0_title);
  strcat(fname12, ".dat");
  print_to_plot(row_l2, fname12); // fname + "/2_b-trans_l2_tca.dat");


  //PAPI Counter2
  char fname13[100];
  strcpy(fname13, fname);
  strcat(fname13, "naive_");
  strcat(fname13, title);
  strcat(fname13, "_");
  strcat(fname13, count1_title);
  strcat(fname13, ".dat");
  print_to_plot(naive_l3, fname13); // fname + "/3_naive_l3_tca.dat");  print_to_plot(naive_l2, fname7); // fname + "/2_naive_l2_tca.dat");

  char fname14[100];
  strcpy(fname14, fname);
  strcat(fname14, "naive2_");
  strcat(fname14, title);
  strcat(fname14, "_");
  strcat(fname14, count1_title);
  strcat(fname14, ".dat");
  print_to_plot(naive2_l3, fname14); // fname + "/3_naive2_l3_tca.dat");

  char fname15[100];
  strcpy(fname15, fname);
  strcat(fname15, "naive2_par2_");
  strcat(fname15, title);
  strcat(fname15, "_");
  strcat(fname15, count1_title);
  strcat(fname15, ".dat");
  print_to_plot(naive2_par2_l3, fname15); // fname + "/3_naive2_par2_l3_tca.dat");

  char fname16[100];
  strcpy(fname16, fname);
  strcat(fname16, "naive2_par4_");
  strcat(fname16, title);
  strcat(fname16, "_");
  strcat(fname16, count1_title);
  strcat(fname16, ".dat");
  print_to_plot(naive2_par4_l3, fname16); // fname + "/3_naive2_par4_l3_tca.dat");

  char fname17[100];
  strcpy(fname17, fname);
  strcat(fname17, "a-trans_");
  strcat(fname17, title);
  strcat(fname17, "_");
  strcat(fname17, count1_title);
  strcat(fname17, ".dat");
  print_to_plot(col_l3, fname17); // fname + "/3_a-trans_l3_tca.dat");

  char fname18[100];
  strcpy(fname18, fname);
  strcat(fname18, "b-trans_");
  strcat(fname18, title);
  strcat(fname18, "_");
  strcat(fname18, count1_title);
  strcat(fname18, ".dat");
  print_to_plot(row_l3, fname18); // "/3_b-trans_l3_trans.dat");
  
}

void test_runtime_square_mirror(int from, int to, int step, int testruns, int events[], int event_size, const char* fname, const char* title, const char* count0_title, const char* count1_title) {

  vector<pair<int, double > > naive_runtime, naive2_runtime, naive2_par2_runtime, naive2_par4_runtime, col_runtime, row_runtime;
  vector<pair<int, double > > naive_l2, naive2_par2_l2, naive2_par4_l2, naive2_l2, col_l2, row_l2;
  vector<pair<int, double > > naive_l3, naive2_par2_l3, naive2_par4_l3, naive2_l3, col_l3, row_l3;

  int size;

  //int testruns = 1;
  //int events[2] = {PAPI_L2_TCA, PAPI_L3_TCA};
  //int event_size = 2;

  for (int run=from; run<to+1; run+=step) {

    size = run;

    cout << "Running test n = " << run << endl;
    
    int m = run;
    int n = run;
    int p = run;
  
    matrix A, B;
  
    resize_matrix(A, m, n);
    resize_matrix(B, n, p);

    cout << " - starting test of row_mult" << endl;
    pair<double,pair<double,double> > res_row = test_row_mult_exclude_preprocess(A, B, events, event_size, testruns);
    row_runtime.push_back(make_pair(run, res_row.first));
    row_l2.push_back(make_pair(run, res_row.second.first));
    row_l3.push_back(make_pair(run, res_row.second.second));

    cout << " - starting test of col_mult" << endl;
    pair<double,pair<double,double> > res_col = test_col_mult_exclude_preprocess(A, B, events, event_size, testruns);
    col_runtime.push_back(make_pair(run, res_col.first));
    col_l2.push_back(make_pair(run, res_col.second.first));
    col_l3.push_back(make_pair(run, res_col.second.second));

    cout << " - starting test of naive_mult" << endl;
    pair<double,pair<double,double> > res_naive = test_naive_mult_exclude_preprocess(A, B, events, event_size, testruns);
    naive_runtime.push_back(make_pair(run, res_naive.first));
    naive_l2.push_back(make_pair(run, res_naive.second.first));
    naive_l3.push_back(make_pair(run, res_naive.second.second));    
    
    cout << " - starting test of naive2_mult" << endl;
    pair<double,pair<double,double> > res_naive2 = test_naive2_mult_exclude_preprocess(size, events, event_size, testruns);
    naive2_runtime.push_back(make_pair(run, res_naive2.first));
    naive2_l2.push_back(make_pair(run, res_naive2.second.first));
    naive2_l3.push_back(make_pair(run, res_naive2.second.second));    

    cout << " - starting test of naive2_par2_mult" << endl;
    pair<double,pair<double,double> > res_naive2_par2 = test_naive2_parallel_mult_exclude_preprocess(size, 2, events, event_size, testruns);
    naive2_par2_runtime.push_back(make_pair(run, res_naive2_par2.first));
    naive2_par2_l2.push_back(make_pair(run, res_naive2_par2.second.first));
    naive2_par2_l3.push_back(make_pair(run, res_naive2_par2.second.second));    

    cout << " - starting test of naive2_par4_mult" << endl;
    pair<double,pair<double,double> > res_naive2_par4 = test_naive2_parallel_mult_exclude_preprocess(size, 4, events, event_size, testruns);
    naive2_par4_runtime.push_back(make_pair(run, res_naive2_par4.first));
    naive2_par4_l2.push_back(make_pair(run, res_naive2_par4.second.first));
    naive2_par4_l3.push_back(make_pair(run, res_naive2_par4.second.second));

  }

  char fname1[100];
  strcpy(fname1,fname);
  strcat(fname1, "naive_");
  strcat(fname1, title);
  strcat(fname1, "_runtime.dat");
  print_to_plot(naive_runtime, fname1); // strcat(fname, "/1_naive_l2_l3_tca_runtime.dat"));

  char fname2[100];
  strcpy(fname2,fname);
  strcat(fname2, "naive2_");
  strcat(fname2, title);
  strcat(fname2, "_runtime.dat");
  print_to_plot(naive2_runtime, fname2); // fname + "/1_naive2_l2_l3_tca_runtime.dat");

  char fname3[100];
  strcpy(fname3,fname);
  strcat(fname3, "naive2_par2_");
  strcat(fname3, title);
  strcat(fname3, "_runtime.dat");
  print_to_plot(naive2_par2_runtime, fname3); // fname + "/1_naive2_par2_l2_l3_tca_runtime.dat");

  char fname4[100];
  strcpy(fname4,fname);
  strcat(fname4, "naive2_par4_");
  strcat(fname4, title);
  strcat(fname4, "_runtime.dat");
  print_to_plot(naive2_par4_runtime, fname4); // fname + "/1_naive2_par4_l2_l3_tca_runtime.dat");

  char fname5[100];
  strcpy(fname5,fname);
  strcat(fname5, "a-trans_");
  strcat(fname5, title);
  strcat(fname5, "_runtime.dat");
  print_to_plot(col_runtime, fname5); // fname + "/1_a-trans_l2_l3_tca_runtime.dat");

  char fname6[100];
  strcpy(fname6,fname);
  strcat(fname6, "b-trans_");
  strcat(fname6, title);
  strcat(fname6, "_runtime.dat");
  print_to_plot(row_runtime, fname6); // fname + "/1_b-trans_l2_l3_tca_runtime.dat");

  // PAPI Counter1
  char fname7[100];
  strcpy(fname7, fname);
  strcat(fname7, "naive_");
  strcat(fname7, title);
  strcat(fname7, "_");
  strcat(fname7, count0_title);
  strcat(fname7, ".dat");
  print_to_plot(naive_l2, fname7); // fname + "/2_naive_l2_tca.dat");

  char fname8[100];
  strcpy(fname8, fname);
  strcat(fname8, "naive2_");
  strcat(fname8, title);
  strcat(fname8, "_");
  strcat(fname8, count0_title);
  strcat(fname8, ".dat");
  print_to_plot(naive2_l2, fname8); // fname + "/2_naive2_l2_tca.dat");

  char fname9[100];
  strcpy(fname9, fname);
  strcat(fname9, "naive2_par2_");
  strcat(fname9, title);
  strcat(fname9, "_");
  strcat(fname9, count0_title);
  strcat(fname9, ".dat");
  print_to_plot(naive2_par2_l2, fname9); // fname + "/2_naive2_par2_l2_tca.dat");

  char fname10[100];
  strcpy(fname10, fname);
  strcat(fname10, "naive2_par4_");
  strcat(fname10, title);
  strcat(fname10, "_");
  strcat(fname10, count0_title);
  strcat(fname10, ".dat");
  print_to_plot(naive2_par4_l2, fname10); // fname + "/2_naive_par4_l2_tca.dat");

  char fname11[100];
  strcpy(fname11, fname);
  strcat(fname11, "a-trans_");
  strcat(fname11, title);
  strcat(fname11, "_");
  strcat(fname11, count0_title);
  strcat(fname11, ".dat");
  print_to_plot(col_l2, fname11); // fname + "/2_a-trans_l2_tca.dat");

  char fname12[100];
  strcpy(fname12, fname);
  strcat(fname12, "b-trans_");
  strcat(fname12, title);
  strcat(fname12, "_");
  strcat(fname12, count0_title);
  strcat(fname12, ".dat");
  print_to_plot(row_l2, fname12); // fname + "/2_b-trans_l2_tca.dat");


  //PAPI Counter2
  char fname13[100];
  strcpy(fname13, fname);
  strcat(fname13, "naive_");
  strcat(fname13, title);
  strcat(fname13, "_");
  strcat(fname13, count1_title);
  strcat(fname13, ".dat");
  print_to_plot(naive_l3, fname13); // fname + "/3_naive_l3_tca.dat");  print_to_plot(naive_l2, fname7); // fname + "/2_naive_l2_tca.dat");

  char fname14[100];
  strcpy(fname14, fname);
  strcat(fname14, "naive2_");
  strcat(fname14, title);
  strcat(fname14, "_");
  strcat(fname14, count1_title);
  strcat(fname14, ".dat");
  print_to_plot(naive2_l3, fname14); // fname + "/3_naive2_l3_tca.dat");

  char fname15[100];
  strcpy(fname15, fname);
  strcat(fname15, "naive2_par2_");
  strcat(fname15, title);
  strcat(fname15, "_");
  strcat(fname15, count1_title);
  strcat(fname15, ".dat");
  print_to_plot(naive2_par2_l3, fname15); // fname + "/3_naive2_par2_l3_tca.dat");

  char fname16[100];
  strcpy(fname16, fname);
  strcat(fname16, "naive2_par4_");
  strcat(fname16, title);
  strcat(fname16, "_");
  strcat(fname16, count1_title);
  strcat(fname16, ".dat");
  print_to_plot(naive2_par4_l3, fname16); // fname + "/3_naive2_par4_l3_tca.dat");

  char fname17[100];
  strcpy(fname17, fname);
  strcat(fname17, "a-trans_");
  strcat(fname17, title);
  strcat(fname17, "_");
  strcat(fname17, count1_title);
  strcat(fname17, ".dat");
  print_to_plot(col_l3, fname17); // fname + "/3_a-trans_l3_tca.dat");

  char fname18[100];
  strcpy(fname18, fname);
  strcat(fname18, "b-trans_");
  strcat(fname18, title);
  strcat(fname18, "_");
  strcat(fname18, count1_title);
  strcat(fname18, ".dat");
  print_to_plot(row_l3, fname18); // "/3_b-trans_l3_trans.dat");
  
}

void test_runtime_ratio(int from, int to, int testruns, int events[], int event_size, const char* fname, const char* title, const char* count_title) {

  vector<pair<int, double > > naive_runtime, naive2_runtime, naive2_par2_runtime, naive2_par4_runtime, col_runtime, row_runtime;
  vector<pair<int, double > > naive_l2, naive2_par2_l2, naive2_par4_l2, naive2_l2, col_l2, row_l2;
  vector<pair<int, double > > naive_l3, naive2_par2_l3, naive2_par4_l3, naive2_l3, col_l3, row_l3;

  int size;

  //int testruns = 1;
  //int events[2] = {PAPI_L2_TCA, PAPI_L3_TCA};
  //int event_size = 2;

  for (int run=from; run<to+1; run++) {

    size = run;

    cout << "Running test n = " << run << endl;
    
    int m = run;
    int n = 4096/4;
    int p = run;
  
    matrix A, B;
  
    resize_matrix(A, m, n);
    resize_matrix(B, n, p);

    cout << " - starting test of row_mult" << endl;
    pair<double,pair<double,double> > res_row = test_row_mult_exclude_preprocess(A, B, events, event_size, testruns);
    row_runtime.push_back(make_pair(run, res_row.first));
    row_l2.push_back(make_pair(run, res_row.second.first/res_row.second.second));
    row_l3.push_back(make_pair(run, res_row.second.second));

    if (run <= 200) {
      cout << " - starting test of col_mult" << endl;
      pair<double,pair<double,double> > res_col = test_col_mult_exclude_preprocess(A, B, events, event_size, testruns);
      col_runtime.push_back(make_pair(run, res_col.first));
      col_l2.push_back(make_pair(run, res_col.second.first/res_col.second.second));
      col_l3.push_back(make_pair(run, res_col.second.second));
    }

    if (run <= 200) {
    cout << " - starting test of naive_mult" << endl;
    pair<double,pair<double,double> > res_naive = test_naive_mult_exclude_preprocess(A, B, events, event_size, testruns);
    naive_runtime.push_back(make_pair(run, res_naive.first));
    naive_l2.push_back(make_pair(run, res_naive.second.first/res_naive.second.second));
    naive_l3.push_back(make_pair(run, res_naive.second.second));    
    }
    
    cout << " - starting test of naive2_mult" << endl;
    pair<double,pair<double,double> > res_naive2 = test_naive2_mult_exclude_preprocess(size, events, event_size, testruns);
    naive2_runtime.push_back(make_pair(run, res_naive2.first));
    naive2_l2.push_back(make_pair(run, res_naive2.second.first/res_naive2.second.second));
    naive2_l3.push_back(make_pair(run, res_naive2.second.second));    

    cout << " - starting test of naive2_par2_mult" << endl;
    pair<double,pair<double,double> > res_naive2_par2 = test_naive2_parallel_mult_exclude_preprocess(size, 2, events, event_size, testruns);
    naive2_par2_runtime.push_back(make_pair(run, res_naive2_par2.first));
    naive2_par2_l2.push_back(make_pair(run, res_naive2_par2.second.first/res_naive2_par2.second.second));
    naive2_par2_l3.push_back(make_pair(run, res_naive2_par2.second.second));    

    cout << " - starting test of naive2_par4_mult" << endl;
    pair<double,pair<double,double> > res_naive2_par4 = test_naive2_parallel_mult_exclude_preprocess(size, 4, events, event_size, testruns);
    naive2_par4_runtime.push_back(make_pair(run, res_naive2_par4.first));
    naive2_par4_l2.push_back(make_pair(run, res_naive2_par4.second.first/res_naive2_par4.second.second));
    naive2_par4_l3.push_back(make_pair(run, res_naive2_par4.second.second));

  }

  char fname1[100];
  strcpy(fname1,fname);
  strcat(fname1, "naive_");
  strcat(fname1, title);
  strcat(fname1, "_runtime.dat");
  print_to_plot(naive_runtime, fname1); // strcat(fname, "/1_naive_l2_l3_tca_runtime.dat"));

  char fname2[100];
  strcpy(fname2,fname);
  strcat(fname2, "naive2_");
  strcat(fname2, title);
  strcat(fname2, "_runtime.dat");
  print_to_plot(naive2_runtime, fname2); // fname + "/1_naive2_l2_l3_tca_runtime.dat");

  char fname3[100];
  strcpy(fname3,fname);
  strcat(fname3, "naive2_par2_");
  strcat(fname3, title);
  strcat(fname3, "_runtime.dat");
  print_to_plot(naive2_par2_runtime, fname3); // fname + "/1_naive2_par2_l2_l3_tca_runtime.dat");

  char fname4[100];
  strcpy(fname4,fname);
  strcat(fname4, "naive2_par4_");
  strcat(fname4, title);
  strcat(fname4, "_runtime.dat");
  print_to_plot(naive2_par4_runtime, fname4); // fname + "/1_naive2_par4_l2_l3_tca_runtime.dat");

  char fname5[100];
  strcpy(fname5,fname);
  strcat(fname5, "a-trans_");
  strcat(fname5, title);
  strcat(fname5, "_runtime.dat");
  print_to_plot(col_runtime, fname5); // fname + "/1_a-trans_l2_l3_tca_runtime.dat");

  char fname6[100];
  strcpy(fname6,fname);
  strcat(fname6, "b-trans_");
  strcat(fname6, title);
  strcat(fname6, "_runtime.dat");
  print_to_plot(row_runtime, fname6); // fname + "/1_b-trans_l2_l3_tca_runtime.dat");

  // PAPI Counter1
  char fname7[100];
  strcpy(fname7, fname);
  strcat(fname7, "naive_");
  strcat(fname7, title);
  strcat(fname7, "_");
  strcat(fname7, count_title);
  strcat(fname7, ".dat");
  print_to_plot(naive_l2, fname7); // fname + "/2_naive_l2_tca.dat");

  char fname8[100];
  strcpy(fname8, fname);
  strcat(fname8, "naive2_");
  strcat(fname8, title);
  strcat(fname8, "_");
  strcat(fname8, count_title);
  strcat(fname8, ".dat");
  print_to_plot(naive2_l2, fname8); // fname + "/2_naive2_l2_tca.dat");

  char fname9[100];
  strcpy(fname9, fname);
  strcat(fname9, "naive2_par2_");
  strcat(fname9, title);
  strcat(fname9, "_");
  strcat(fname9, count_title);
  strcat(fname9, ".dat");
  print_to_plot(naive2_par2_l2, fname9); // fname + "/2_naive2_par2_l2_tca.dat");

  char fname10[100];
  strcpy(fname10, fname);
  strcat(fname10, "naive2_par4_");
  strcat(fname10, title);
  strcat(fname10, "_");
  strcat(fname10, count_title);
  strcat(fname10, ".dat");
  print_to_plot(naive2_par4_l2, fname10); // fname + "/2_naive_par4_l2_tca.dat");

  char fname11[100];
  strcpy(fname11, fname);
  strcat(fname11, "a-trans_");
  strcat(fname11, title);
  strcat(fname11, "_");
  strcat(fname11, count_title);
  strcat(fname11, ".dat");
  print_to_plot(col_l2, fname11); // fname + "/2_a-trans_l2_tca.dat");

  char fname12[100];
  strcpy(fname12, fname);
  strcat(fname12, "b-trans_");
  strcat(fname12, title);
  strcat(fname12, "_");
  strcat(fname12, count_title);
  strcat(fname12, ".dat");
  print_to_plot(row_l2, fname12); // fname + "/2_b-trans_l2_tca.dat");


/*
  //PAPI Counter2
  char fname13[100];
  strcpy(fname13, fname);
  strcat(fname13, "naive_");
  strcat(fname13, title);
  strcat(fname13, "_");
  strcat(fname13, count1_title);
  strcat(fname13, ".dat");
  print_to_plot(naive_l3, fname13); // fname + "/3_naive_l3_tca.dat");  print_to_plot(naive_l2, fname7); // fname + "/2_naive_l2_tca.dat");

  char fname14[100];
  strcpy(fname14, fname);
  strcat(fname14, "naive2_");
  strcat(fname14, title);
  strcat(fname14, "_");
  strcat(fname14, count1_title);
  strcat(fname14, ".dat");
  print_to_plot(naive2_l3, fname14); // fname + "/3_naive2_l3_tca.dat");

  char fname15[100];
  strcpy(fname15, fname);
  strcat(fname15, "naive2_par2_");
  strcat(fname15, title);
  strcat(fname15, "_");
  strcat(fname15, count1_title);
  strcat(fname15, ".dat");
  print_to_plot(naive2_par2_l3, fname15); // fname + "/3_naive2_par2_l3_tca.dat");

  char fname16[100];
  strcpy(fname16, fname);
  strcat(fname16, "naive2_par4_");
  strcat(fname16, title);
  strcat(fname16, "_");
  strcat(fname16, count1_title);
  strcat(fname16, ".dat");
  print_to_plot(naive2_par4_l3, fname16); // fname + "/3_naive2_par4_l3_tca.dat");

  char fname17[100];
  strcpy(fname17, fname);
  strcat(fname17, "a-trans_");
  strcat(fname17, title);
  strcat(fname17, "_");
  strcat(fname17, count1_title);
  strcat(fname17, ".dat");
  print_to_plot(col_l3, fname17); // fname + "/3_a-trans_l3_tca.dat");

  char fname18[100];
  strcpy(fname18, fname);
  strcat(fname18, "b-trans_");
  strcat(fname18, title);
  strcat(fname18, "_");
  strcat(fname18, count1_title);
  strcat(fname18, ".dat");
  print_to_plot(row_l3, fname18); // "/3_b-trans_l3_trans.dat");
*/

}

void test_runtime_ratio_square(int from, int to, int step, int testruns, int events[], int event_size, const char* fname, const char* title, const char* count_title) {

  vector<pair<int, double > > naive_runtime, naive2_runtime, naive2_par2_runtime, naive2_par4_runtime, col_runtime, row_runtime;
  vector<pair<int, double > > naive_l2, naive2_par2_l2, naive2_par4_l2, naive2_l2, col_l2, row_l2;
  vector<pair<int, double > > naive_l3, naive2_par2_l3, naive2_par4_l3, naive2_l3, col_l3, row_l3;

  int size;

  //int testruns = 1;
  //int events[2] = {PAPI_L2_TCA, PAPI_L3_TCA};
  //int event_size = 2;

  for (int run=from; run<to+1; run+=step) {

    size = run;

    cout << "Running test n = " << run << endl;
    
    int m = run;
    int n = run;
    int p = run;
  
    matrix A, B;
  
    resize_matrix(A, m, n);
    resize_matrix(B, n, p);

    cout << " - starting test of row_mult" << endl;
    pair<double,pair<double,double> > res_row = test_row_mult_exclude_preprocess(A, B, events, event_size, testruns);
    row_runtime.push_back(make_pair(run, res_row.first));
    row_l2.push_back(make_pair(run, res_row.second.first/res_row.second.second));
    row_l3.push_back(make_pair(run, res_row.second.second));

    if (run <= 200) {
      cout << " - starting test of col_mult" << endl;
      pair<double,pair<double,double> > res_col = test_col_mult_exclude_preprocess(A, B, events, event_size, testruns);
      col_runtime.push_back(make_pair(run, res_col.first));
      col_l2.push_back(make_pair(run, res_col.second.first/res_col.second.second));
      col_l3.push_back(make_pair(run, res_col.second.second));
    }

    if (run <= 200) {
    cout << " - starting test of naive_mult" << endl;
    pair<double,pair<double,double> > res_naive = test_naive_mult_exclude_preprocess(A, B, events, event_size, testruns);
    naive_runtime.push_back(make_pair(run, res_naive.first));
    naive_l2.push_back(make_pair(run, res_naive.second.first/res_naive.second.second));
    naive_l3.push_back(make_pair(run, res_naive.second.second));    
    }
    
    cout << " - starting test of naive2_mult" << endl;
    pair<double,pair<double,double> > res_naive2 = test_naive2_mult_exclude_preprocess(size, events, event_size, testruns);
    naive2_runtime.push_back(make_pair(run, res_naive2.first));
    naive2_l2.push_back(make_pair(run, res_naive2.second.first/res_naive2.second.second));
    naive2_l3.push_back(make_pair(run, res_naive2.second.second));    

    cout << " - starting test of naive2_par2_mult" << endl;
    pair<double,pair<double,double> > res_naive2_par2 = test_naive2_parallel_mult_exclude_preprocess(size, 2, events, event_size, testruns);
    naive2_par2_runtime.push_back(make_pair(run, res_naive2_par2.first));
    naive2_par2_l2.push_back(make_pair(run, res_naive2_par2.second.first/res_naive2_par2.second.second));
    naive2_par2_l3.push_back(make_pair(run, res_naive2_par2.second.second));    

    cout << " - starting test of naive2_par4_mult" << endl;
    pair<double,pair<double,double> > res_naive2_par4 = test_naive2_parallel_mult_exclude_preprocess(size, 4, events, event_size, testruns);
    naive2_par4_runtime.push_back(make_pair(run, res_naive2_par4.first));
    naive2_par4_l2.push_back(make_pair(run, res_naive2_par4.second.first/res_naive2_par4.second.second));
    naive2_par4_l3.push_back(make_pair(run, res_naive2_par4.second.second));

  }

  char fname1[100];
  strcpy(fname1,fname);
  strcat(fname1, "naive_");
  strcat(fname1, title);
  strcat(fname1, "_runtime.dat");
  print_to_plot(naive_runtime, fname1); // strcat(fname, "/1_naive_l2_l3_tca_runtime.dat"));

  char fname2[100];
  strcpy(fname2,fname);
  strcat(fname2, "naive2_");
  strcat(fname2, title);
  strcat(fname2, "_runtime.dat");
  print_to_plot(naive2_runtime, fname2); // fname + "/1_naive2_l2_l3_tca_runtime.dat");

  char fname3[100];
  strcpy(fname3,fname);
  strcat(fname3, "naive2_par2_");
  strcat(fname3, title);
  strcat(fname3, "_runtime.dat");
  print_to_plot(naive2_par2_runtime, fname3); // fname + "/1_naive2_par2_l2_l3_tca_runtime.dat");

  char fname4[100];
  strcpy(fname4,fname);
  strcat(fname4, "naive2_par4_");
  strcat(fname4, title);
  strcat(fname4, "_runtime.dat");
  print_to_plot(naive2_par4_runtime, fname4); // fname + "/1_naive2_par4_l2_l3_tca_runtime.dat");

  char fname5[100];
  strcpy(fname5,fname);
  strcat(fname5, "a-trans_");
  strcat(fname5, title);
  strcat(fname5, "_runtime.dat");
  print_to_plot(col_runtime, fname5); // fname + "/1_a-trans_l2_l3_tca_runtime.dat");

  char fname6[100];
  strcpy(fname6,fname);
  strcat(fname6, "b-trans_");
  strcat(fname6, title);
  strcat(fname6, "_runtime.dat");
  print_to_plot(row_runtime, fname6); // fname + "/1_b-trans_l2_l3_tca_runtime.dat");

  // PAPI Counter1
  char fname7[100];
  strcpy(fname7, fname);
  strcat(fname7, "naive_");
  strcat(fname7, title);
  strcat(fname7, "_");
  strcat(fname7, count_title);
  strcat(fname7, ".dat");
  print_to_plot(naive_l2, fname7); // fname + "/2_naive_l2_tca.dat");

  char fname8[100];
  strcpy(fname8, fname);
  strcat(fname8, "naive2_");
  strcat(fname8, title);
  strcat(fname8, "_");
  strcat(fname8, count_title);
  strcat(fname8, ".dat");
  print_to_plot(naive2_l2, fname8); // fname + "/2_naive2_l2_tca.dat");

  char fname9[100];
  strcpy(fname9, fname);
  strcat(fname9, "naive2_par2_");
  strcat(fname9, title);
  strcat(fname9, "_");
  strcat(fname9, count_title);
  strcat(fname9, ".dat");
  print_to_plot(naive2_par2_l2, fname9); // fname + "/2_naive2_par2_l2_tca.dat");

  char fname10[100];
  strcpy(fname10, fname);
  strcat(fname10, "naive2_par4_");
  strcat(fname10, title);
  strcat(fname10, "_");
  strcat(fname10, count_title);
  strcat(fname10, ".dat");
  print_to_plot(naive2_par4_l2, fname10); // fname + "/2_naive_par4_l2_tca.dat");

  char fname11[100];
  strcpy(fname11, fname);
  strcat(fname11, "a-trans_");
  strcat(fname11, title);
  strcat(fname11, "_");
  strcat(fname11, count_title);
  strcat(fname11, ".dat");
  print_to_plot(col_l2, fname11); // fname + "/2_a-trans_l2_tca.dat");

  char fname12[100];
  strcpy(fname12, fname);
  strcat(fname12, "b-trans_");
  strcat(fname12, title);
  strcat(fname12, "_");
  strcat(fname12, count_title);
  strcat(fname12, ".dat");
  print_to_plot(row_l2, fname12); // fname + "/2_b-trans_l2_tca.dat");


/*
  //PAPI Counter2
  char fname13[100];
  strcpy(fname13, fname);
  strcat(fname13, "naive_");
  strcat(fname13, title);
  strcat(fname13, "_");
  strcat(fname13, count1_title);
  strcat(fname13, ".dat");
  print_to_plot(naive_l3, fname13); // fname + "/3_naive_l3_tca.dat");  print_to_plot(naive_l2, fname7); // fname + "/2_naive_l2_tca.dat");

  char fname14[100];
  strcpy(fname14, fname);
  strcat(fname14, "naive2_");
  strcat(fname14, title);
  strcat(fname14, "_");
  strcat(fname14, count1_title);
  strcat(fname14, ".dat");
  print_to_plot(naive2_l3, fname14); // fname + "/3_naive2_l3_tca.dat");

  char fname15[100];
  strcpy(fname15, fname);
  strcat(fname15, "naive2_par2_");
  strcat(fname15, title);
  strcat(fname15, "_");
  strcat(fname15, count1_title);
  strcat(fname15, ".dat");
  print_to_plot(naive2_par2_l3, fname15); // fname + "/3_naive2_par2_l3_tca.dat");

  char fname16[100];
  strcpy(fname16, fname);
  strcat(fname16, "naive2_par4_");
  strcat(fname16, title);
  strcat(fname16, "_");
  strcat(fname16, count1_title);
  strcat(fname16, ".dat");
  print_to_plot(naive2_par4_l3, fname16); // fname + "/3_naive2_par4_l3_tca.dat");

  char fname17[100];
  strcpy(fname17, fname);
  strcat(fname17, "a-trans_");
  strcat(fname17, title);
  strcat(fname17, "_");
  strcat(fname17, count1_title);
  strcat(fname17, ".dat");
  print_to_plot(col_l3, fname17); // fname + "/3_a-trans_l3_tca.dat");

  char fname18[100];
  strcpy(fname18, fname);
  strcat(fname18, "b-trans_");
  strcat(fname18, title);
  strcat(fname18, "_");
  strcat(fname18, count1_title);
  strcat(fname18, ".dat");
  print_to_plot(row_l3, fname18); // "/3_b-trans_l3_trans.dat");
*/

}

void test_runtime_l2_l3_accesses_exp2() {

  vector<pair<int, double > > naive_runtime, naive2_runtime, naive2_par2_runtime, naive2_par4_runtime, col_runtime, row_runtime;
  vector<pair<int, double > > naive_l2, naive2_par2_l2, naive2_par4_l2, naive2_l2, col_l2, row_l2;
  vector<pair<int, double > > naive_l3, naive2_par2_l3, naive2_par4_l3, naive2_l3, col_l3, row_l3;

  int size;

  int testruns = 1;
  int events[2] = {PAPI_TLB_DM, PAPI_TOT_CYC};
  int event_size = 2;

  for (int run=2; run<501; run++) {

    size = run;

    cout << "Running test n = " << run << endl;
    
    int m = run;
    int n = 4096;
    int p = run;
  
    matrix A, B;
  
    resize_matrix(A, m, n);
    resize_matrix(B, n, p);

    cout << " - starting test of row_mult" << endl;
    pair<double,pair<double,double> > res_row = test_row_mult_exclude_preprocess(A, B, events, event_size, testruns);
    row_runtime.push_back(make_pair(run, res_row.first));
    row_l2.push_back(make_pair(run, res_row.second.first));
    row_l3.push_back(make_pair(run, res_row.second.second));

    if (run <= 200) {
      cout << " - starting test of col_mult" << endl;
      pair<double,pair<double,double> > res_col = test_col_mult_exclude_preprocess(A, B, events, event_size, testruns);
      col_runtime.push_back(make_pair(run, res_col.first));
      col_l2.push_back(make_pair(run, res_col.second.first));
      col_l3.push_back(make_pair(run, res_col.second.second));
    }

    if (run <= 200) {
    cout << " - starting test of naive_mult" << endl;
    pair<double,pair<double,double> > res_naive = test_naive_mult_exclude_preprocess(A, B, events, event_size, testruns);
    naive_runtime.push_back(make_pair(run, res_naive.first));
    naive_l2.push_back(make_pair(run, res_naive.second.first));
    naive_l3.push_back(make_pair(run, res_naive.second.second));    
    }
    
    cout << " - starting test of naive2_mult" << endl;
    pair<double,pair<double,double> > res_naive2 = test_naive2_mult_exclude_preprocess(size, events, event_size, testruns);
    naive2_runtime.push_back(make_pair(run, res_naive2.first));
    naive2_l2.push_back(make_pair(run, res_naive2.second.first));
    naive2_l3.push_back(make_pair(run, res_naive2.second.second));    

    cout << " - starting test of naive2_par2_mult" << endl;
    pair<double,pair<double,double> > res_naive2_par2 = test_naive2_parallel_mult_exclude_preprocess(size, 2, events, event_size, testruns);
    naive2_par2_runtime.push_back(make_pair(run, res_naive2_par2.first));
    naive2_par2_l2.push_back(make_pair(run, res_naive2_par2.second.first));
    naive2_par2_l3.push_back(make_pair(run, res_naive2_par2.second.second));    

    cout << " - starting test of naive2_par4_mult" << endl;
    pair<double,pair<double,double> > res_naive2_par4 = test_naive2_parallel_mult_exclude_preprocess(size, 4, events, event_size, testruns);
    naive2_par4_runtime.push_back(make_pair(run, res_naive2_par4.first));
    naive2_par4_l2.push_back(make_pair(run, res_naive2_par4.second.first));
    naive2_par4_l3.push_back(make_pair(run, res_naive2_par4.second.second));

  }

  print_to_plot(naive_runtime, "data_matrix/Experiment2/naive_runtime.dat");
  print_to_plot(naive2_runtime, "data_matrix/Experiment2/naive2_runtime.dat");
  print_to_plot(naive2_par2_runtime, "data_matrix/Experiment2/naive2_par2_runtime.dat");
  print_to_plot(naive2_par4_runtime, "data_matrix/Experiment2/naive2_par4_runtime.dat");
  print_to_plot(col_runtime, "data_matrix/Experiment2/col_runtime.dat");
  print_to_plot(row_runtime, "data_matrix/Experiment2/row_runtime.dat");

  print_to_plot(naive_l2, "data_matrix/Experiment2/naive_tlb.dat");
  print_to_plot(naive2_l2, "data_matrix/Experiment2/naive2_tlb.dat");
  print_to_plot(naive2_par2_l2, "data_matrix/Experiment2/naive2_par2_tlb.dat");
  print_to_plot(naive2_par4_l2, "data_matrix/Experiment2/naive_par4_tlb.dat");
  print_to_plot(col_l2, "data_matrix/Experiment2/col_tlb.dat");
  print_to_plot(row_l2, "data_matrix/Experiment2/row_tlb.dat");

  print_to_plot(naive_l3, "data_matrix/Experiment2/naive_cyc.dat");
  print_to_plot(naive2_l3, "data_matrix/Experiment2/naive2_cyc.dat");
  print_to_plot(naive2_par2_l3, "data_matrix/Experiment2/naive_par2_cyc.dat");
  print_to_plot(naive2_par4_l3, "data_matrix/Experiment2/naive_par4_cyc.dat");
  print_to_plot(col_l3, "data_matrix/Experiment2/col_cyc.dat");
  print_to_plot(row_l3, "data_matrix/Experiment2/row_cyc.dat");

}


void test_runtime_l2_l3_accesses_exp3() {

  vector<pair<int, double > > naive_runtime, naive2_runtime, naive2_par2_runtime, naive2_par4_runtime, col_runtime, row_runtime;
  vector<pair<int, double > > naive_l2, naive2_par2_l2, naive2_par4_l2, naive2_l2, col_l2, row_l2;
  vector<pair<int, double > > naive_l3, naive2_par2_l3, naive2_par4_l3, naive2_l3, col_l3, row_l3;

  int size;

  int testruns = 1;
  int events[2] = {PAPI_L2_TCA, PAPI_L2_TCM};
  int event_size = 2;

  for (int run=2; run<501; run++) {

    size = run;

    cout << "Running test n = " << run << endl;
    
    int m = run;
    int n = 4096;
    int p = run;
  
    matrix A, B;
  
    resize_matrix(A, m, n);
    resize_matrix(B, n, p);

    cout << " - starting test of row_mult" << endl;
    pair<double,pair<double,double> > res_row = test_row_mult_exclude_preprocess(A, B, events, event_size, testruns);
    row_runtime.push_back(make_pair(run, res_row.first));
    row_l2.push_back(make_pair(run, res_row.second.first/res_row.second.second));
    row_l3.push_back(make_pair(run, res_row.second.second));

    if (run <= 200) {
      cout << " - starting test of col_mult" << endl;
      pair<double,pair<double,double> > res_col = test_col_mult_exclude_preprocess(A, B, events, event_size, testruns);
      col_runtime.push_back(make_pair(run, res_col.first));
      col_l2.push_back(make_pair(run, res_col.second.first/res_col.second.second));
      col_l3.push_back(make_pair(run, res_col.second.second));
    }

    if (run <= 200) {
    cout << " - starting test of naive_mult" << endl;
    pair<double,pair<double,double> > res_naive = test_naive_mult_exclude_preprocess(A, B, events, event_size, testruns);
    naive_runtime.push_back(make_pair(run, res_naive.first));
    naive_l2.push_back(make_pair(run, res_naive.second.first/res_naive.second.second));
    naive_l3.push_back(make_pair(run, res_naive.second.second));    
    }
    
    cout << " - starting test of naive2_mult" << endl;
    pair<double,pair<double,double> > res_naive2 = test_naive2_mult_exclude_preprocess(size, events, event_size, testruns);
    naive2_runtime.push_back(make_pair(run, res_naive2.first));
    naive2_l2.push_back(make_pair(run, res_naive2.second.first/res_naive2.second.second));
    naive2_l3.push_back(make_pair(run, res_naive2.second.second));    

    cout << " - starting test of naive2_par2_mult" << endl;
    pair<double,pair<double,double> > res_naive2_par2 = test_naive2_parallel_mult_exclude_preprocess(size, 2, events, event_size, testruns);
    naive2_par2_runtime.push_back(make_pair(run, res_naive2_par2.first));
    naive2_par2_l2.push_back(make_pair(run, res_naive2_par2.second.first/res_naive2_par2.second.second));
    naive2_par2_l3.push_back(make_pair(run, res_naive2_par2.second.second));    

    cout << " - starting test of naive2_par4_mult" << endl;
    pair<double,pair<double,double> > res_naive2_par4 = test_naive2_parallel_mult_exclude_preprocess(size, 4, events, event_size, testruns);
    naive2_par4_runtime.push_back(make_pair(run, res_naive2_par4.first));
    naive2_par4_l2.push_back(make_pair(run, res_naive2_par4.second.first/res_naive2_par4.second.second));
    naive2_par4_l3.push_back(make_pair(run, res_naive2_par4.second.second));

  }

  print_to_plot(naive_runtime, "data_matrix/Experiment3/naive_runtime.dat");
  print_to_plot(naive2_runtime, "data_matrix/Experiment3/naive2_runtime.dat");
  print_to_plot(naive2_par2_runtime, "data_matrix/Experiment3/naive2_par2_runtime.dat");
  print_to_plot(naive2_par4_runtime, "data_matrix/Experiment3/naive2_par4_runtime.dat");
  print_to_plot(col_runtime, "data_matrix/Experiment3/col_runtime.dat");
  print_to_plot(row_runtime, "data_matrix/Experiment3/row_runtime.dat");

  //print_to_plot(naive_l2, "data_matrix/Experiment3/naive_tlb.dat");
  //print_to_plot(naive2_l2, "data_matrix/Experiment3/naive2_tlb.dat");
  //print_to_plot(naive2_par2_l2, "data_matrix/Experiment3/naive2_par2_tlb.dat");
  //print_to_plot(naive2_par4_l2, "data_matrix/Experiment3/naive_par4_tlb.dat");
  //print_to_plot(col_l2, "data_matrix/Experiment3/col_tlb.dat");
  //print_to_plot(row_l2, "data_matrix/Experiment3/row_tlb.dat");

  print_to_plot(naive_l3, "data_matrix/Experiment3/naive_l2_mr.dat");
  print_to_plot(naive2_l3, "data_matrix/Experiment3/naive2_l2_mr.dat");
  print_to_plot(naive2_par2_l3, "data_matrix/Experiment3/naive_par2_mr.dat");
  print_to_plot(naive2_par4_l3, "data_matrix/Experiment3/naive_par4_mr.dat");
  print_to_plot(col_l3, "data_matrix/Experiment3/col_mr.dat");
  print_to_plot(row_l3, "data_matrix/Experiment3/row_mr.dat");

}

void test_runtime_l2_l3_accesses_exp4() {

  vector<pair<int, double > > naive_runtime, naive2_runtime, naive2_par2_runtime, naive2_par4_runtime, col_runtime, row_runtime;
  vector<pair<int, double > > naive_l2, naive2_par2_l2, naive2_par4_l2, naive2_l2, col_l2, row_l2;
  vector<pair<int, double > > naive_l3, naive2_par2_l3, naive2_par4_l3, naive2_l3, col_l3, row_l3;

  int size;

  int testruns = 1;
  int events[2] = {PAPI_TLB_DM, PAPI_TOT_CYC};
  int event_size = 2;

  for (int run=2; run<501; run++) {

    size = run;

    cout << "Running test n = " << run << endl;
    
    int m = run;
    int n = 4096;
    int p = run;
  
    matrix A, B;
  
    resize_matrix(A, m, n);
    resize_matrix(B, n, p);

    cout << " - starting test of row_mult" << endl;
    pair<double,pair<double,double> > res_row = test_row_mult_idx_opt_exclude_preprocess(A, B, events, event_size, testruns);
    row_runtime.push_back(make_pair(run, res_row.first));
    row_l2.push_back(make_pair(run, res_row.second.first));
    row_l3.push_back(make_pair(run, res_row.second.second));

    if (run <= 200) {
      cout << " - starting test of col_mult" << endl;
      pair<double,pair<double,double> > res_col = test_col_mult_exclude_preprocess(A, B, events, event_size, testruns);
      col_runtime.push_back(make_pair(run, res_col.first));
      col_l2.push_back(make_pair(run, res_col.second.first));
      col_l3.push_back(make_pair(run, res_col.second.second));
    }

    if (run <= 200) {
    cout << " - starting test of naive_mult" << endl;
    pair<double,pair<double,double> > res_naive = test_naive_mult_exclude_preprocess(A, B, events, event_size, testruns);
    naive_runtime.push_back(make_pair(run, res_naive.first));
    naive_l2.push_back(make_pair(run, res_naive.second.first));
    naive_l3.push_back(make_pair(run, res_naive.second.second));    
    }
    
    cout << " - starting test of naive2_mult" << endl;
    pair<double,pair<double,double> > res_naive2 = test_naive2_mult_exclude_preprocess(size, events, event_size, testruns);
    naive2_runtime.push_back(make_pair(run, res_naive2.first));
    naive2_l2.push_back(make_pair(run, res_naive2.second.first));
    naive2_l3.push_back(make_pair(run, res_naive2.second.second));    

    cout << " - starting test of naive2_par2_mult" << endl;
    pair<double,pair<double,double> > res_naive2_par2 = test_naive2_parallel_mult_exclude_preprocess(size, 2, events, event_size, testruns);
    naive2_par2_runtime.push_back(make_pair(run, res_naive2_par2.first));
    naive2_par2_l2.push_back(make_pair(run, res_naive2_par2.second.first));
    naive2_par2_l3.push_back(make_pair(run, res_naive2_par2.second.second));    

    cout << " - starting test of naive2_par4_mult" << endl;
    pair<double,pair<double,double> > res_naive2_par4 = test_naive2_parallel_mult_exclude_preprocess(size, 4, events, event_size, testruns);
    naive2_par4_runtime.push_back(make_pair(run, res_naive2_par4.first));
    naive2_par4_l2.push_back(make_pair(run, res_naive2_par4.second.first));
    naive2_par4_l3.push_back(make_pair(run, res_naive2_par4.second.second));

  }

  print_to_plot(naive_runtime, "data_matrix/Experiment4/naive_runtime.dat");
  print_to_plot(naive2_runtime, "data_matrix/Experiment4/naive2_runtime.dat");
  print_to_plot(naive2_par2_runtime, "data_matrix/Experiment4/naive2_par2_runtime.dat");
  print_to_plot(naive2_par4_runtime, "data_matrix/Experiment4/naive2_par4_runtime.dat");
  print_to_plot(col_runtime, "data_matrix/Experiment4/col_runtime.dat");
  print_to_plot(row_runtime, "data_matrix/Experiment4/row_runtime.dat");

  print_to_plot(naive_l2, "data_matrix/Experiment4/naive_tlb.dat");
  print_to_plot(naive2_l2, "data_matrix/Experiment4/naive2_tlb.dat");
  print_to_plot(naive2_par2_l2, "data_matrix/Experiment4/naive2_par2_tlb.dat");
  print_to_plot(naive2_par4_l2, "data_matrix/Experiment4/naive_par4_tlb.dat");
  print_to_plot(col_l2, "data_matrix/Experiment4/col_tlb.dat");
  print_to_plot(row_l2, "data_matrix/Experiment4/row_tlb.dat");

  print_to_plot(naive_l3, "data_matrix/Experiment4/naive_cyc.dat");
  print_to_plot(naive2_l3, "data_matrix/Experiment4/naive2_cyc.dat");
  print_to_plot(naive2_par2_l3, "data_matrix/Experiment4/naive_par2_cyc.dat");
  print_to_plot(naive2_par4_l3, "data_matrix/Experiment4/naive_par4_cyc.dat");
  print_to_plot(col_l3, "data_matrix/Experiment4/col_cyc.dat");
  print_to_plot(row_l3, "data_matrix/Experiment4/row_cyc.dat");

}

void test_runtime_l2_l3_accesses_exp5() {

  vector<pair<int, double > > naive_runtime, naive2_runtime, naive2_par2_runtime, naive2_par4_runtime, col_runtime, row_runtime;
  vector<pair<int, double > > naive_l2, naive2_par2_l2, naive2_par4_l2, naive2_l2, col_l2, row_l2;
  vector<pair<int, double > > naive_l3, naive2_par2_l3, naive2_par4_l3, naive2_l3, col_l3, row_l3;

  int size;

  int testruns = 1;
  int events[2] = {PAPI_L3_TCA, PAPI_L3_TCM};
  int event_size = 2;

  for (int run=2; run<501; run++) {

    size = run;

    cout << "Running test n = " << run << endl;
    
    int m = run;
    int n = 4096;
    int p = run;
  
    matrix A, B;
  
    resize_matrix(A, m, n);
    resize_matrix(B, n, p);

    cout << " - starting test of row_mult" << endl;
    pair<double,pair<double,double> > res_row = test_row_mult_exclude_preprocess(A, B, events, event_size, testruns);
    row_runtime.push_back(make_pair(run, res_row.first));
    row_l2.push_back(make_pair(run, res_row.second.first/res_row.second.second));
    row_l3.push_back(make_pair(run, res_row.second.second));

    if (run <= 200) {
      cout << " - starting test of col_mult" << endl;
      pair<double,pair<double,double> > res_col = test_col_mult_exclude_preprocess(A, B, events, event_size, testruns);
      col_runtime.push_back(make_pair(run, res_col.first));
      col_l2.push_back(make_pair(run, res_col.second.first/res_col.second.second));
      col_l3.push_back(make_pair(run, res_col.second.second));
    }

    if (run <= 200) {
    cout << " - starting test of naive_mult" << endl;
    pair<double,pair<double,double> > res_naive = test_naive_mult_exclude_preprocess(A, B, events, event_size, testruns);
    naive_runtime.push_back(make_pair(run, res_naive.first));
    naive_l2.push_back(make_pair(run, res_naive.second.first/res_naive.second.second));
    naive_l3.push_back(make_pair(run, res_naive.second.second));    
    }
    
    cout << " - starting test of naive2_mult" << endl;
    pair<double,pair<double,double> > res_naive2 = test_naive2_mult_exclude_preprocess(size, events, event_size, testruns);
    naive2_runtime.push_back(make_pair(run, res_naive2.first));
    naive2_l2.push_back(make_pair(run, res_naive2.second.first/res_naive2.second.second));
    naive2_l3.push_back(make_pair(run, res_naive2.second.second));    

    cout << " - starting test of naive2_par2_mult" << endl;
    pair<double,pair<double,double> > res_naive2_par2 = test_naive2_parallel_mult_exclude_preprocess(size, 2, events, event_size, testruns);
    naive2_par2_runtime.push_back(make_pair(run, res_naive2_par2.first));
    naive2_par2_l2.push_back(make_pair(run, res_naive2_par2.second.first/res_naive2_par2.second.second));
    naive2_par2_l3.push_back(make_pair(run, res_naive2_par2.second.second));    

    cout << " - starting test of naive2_par4_mult" << endl;
    pair<double,pair<double,double> > res_naive2_par4 = test_naive2_parallel_mult_exclude_preprocess(size, 4, events, event_size, testruns);
    naive2_par4_runtime.push_back(make_pair(run, res_naive2_par4.first));
    naive2_par4_l2.push_back(make_pair(run, res_naive2_par4.second.first/res_naive2_par4.second.second));
    naive2_par4_l3.push_back(make_pair(run, res_naive2_par4.second.second));

  }

  print_to_plot(naive_runtime, "data_matrix/Experiment5/naive_runtime.dat");
  print_to_plot(naive2_runtime, "data_matrix/Experiment5/naive2_runtime.dat");
  print_to_plot(naive2_par2_runtime, "data_matrix/Experiment5/naive2_par2_runtime.dat");
  print_to_plot(naive2_par4_runtime, "data_matrix/Experiment5/naive2_par4_runtime.dat");
  print_to_plot(col_runtime, "data_matrix/Experiment5/col_runtime.dat");
  print_to_plot(row_runtime, "data_matrix/Experiment5/row_runtime.dat");

  //print_to_plot(naive_l2, "data_matrix/Experiment3/naive_tlb.dat");
  //print_to_plot(naive2_l2, "data_matrix/Experiment3/naive2_tlb.dat");
  //print_to_plot(naive2_par2_l2, "data_matrix/Experiment3/naive2_par2_tlb.dat");
  //print_to_plot(naive2_par4_l2, "data_matrix/Experiment3/naive_par4_tlb.dat");
  //print_to_plot(col_l2, "data_matrix/Experiment3/col_tlb.dat");
  //print_to_plot(row_l2, "data_matrix/Experiment3/row_tlb.dat");

  print_to_plot(naive_l3, "data_matrix/Experiment5/naive_l5_mr.dat");
  print_to_plot(naive2_l3, "data_matrix/Experiment5/naive2_l5_mr.dat");
  print_to_plot(naive2_par2_l3, "data_matrix/Experiment5/naive_par2_mr.dat");
  print_to_plot(naive2_par4_l3, "data_matrix/Experiment5/naive_par4_mr.dat");
  print_to_plot(col_l3, "data_matrix/Experiment5/col_mr.dat");
  print_to_plot(row_l3, "data_matrix/Experiment5/row_mr.dat");

}

void test_runtime_l2_l3_accesses_row() {

  vector<pair<int, double > > naive_runtime, naive2_runtime, naive2_par2_runtime, naive2_par4_runtime, col_runtime, row_runtime;
  vector<pair<int, double > > naive_l2, naive2_par2_l2, naive2_par4_l2, naive2_l2, col_l2, row_l2;
  vector<pair<int, double > > naive_l3, naive2_par2_l3, naive2_par4_l3, naive2_l3, col_l3, row_l3;

  int size;

  int testruns = 1;
  int events[2] = {PAPI_L2_TCA, PAPI_L3_TCA};
  int event_size = 2;

  for (int run=1024; run<4200; run+=256) {

    size = run;

    cout << "Running test on n = " << run << endl;
    
    int m = size;
    int n = size;
    int p = size;
  
    matrix A, B;
  
    resize_matrix(A, m, n);
    resize_matrix(B, n, p);

    cout << " - starting test of row_mult" << endl;
    pair<double,pair<double,double> > res_row = test_row_mult_exclude_preprocess(A, B, events, event_size, testruns);
    row_runtime.push_back(make_pair(run, res_row.first));
    row_l2.push_back(make_pair(run, res_row.second.first));
    row_l3.push_back(make_pair(run, res_row.second.second));

  }

  print_to_plot(row_runtime, "data_matrix/row_runtime_2.dat");
  print_to_plot(row_l2, "data_matrix/row_l2_2.dat");
  print_to_plot(row_l3, "data_matrix/row_l3_2.dat");

}

void test_runtime_l2_l3_accesses_row_rec() {

  vector<pair<int, double > > rec_runtime, row_runtime;
  vector<pair<int, double > > rec_l2, row_l2;
  vector<pair<int, double > > rec_l3, row_l3;

  int size;

  int testruns = 10;
  int numruns = 7;
  int events[2] = {PAPI_L2_TCA, PAPI_L3_TCA};
  int event_size = 2;

  for (int run=1; run<numruns+1; run++) {

    size = 1 << (run-1);

    cout << "Running test " << run << " of " << numruns << endl;
    
    int m = size;
    int n = size;
    int p = size;
  
    matrix A, B;
  
    resize_matrix(A, m, n);
    resize_matrix(B, n, p);

    pair<double,pair<double,double> > res_row = test_row_mult_exclude_preprocess(A, B, events, event_size, testruns);
    row_runtime.push_back(make_pair(size, res_row.first/1000));
    row_l2.push_back(make_pair(size, res_row.second.first));
    row_l3.push_back(make_pair(size, res_row.second.second));

    pair<double,pair<double,double> > res_rec = test_mult_rec(A, B, events, event_size, testruns);
    rec_runtime.push_back(make_pair(size, res_rec.first/1000));
    rec_l2.push_back(make_pair(size, res_rec.second.first));
    rec_l3.push_back(make_pair(size, res_rec.second.second));
   
    //cout << res_naive.first/1000 << " " << (long long) res_naive.second.first << " " << (long long) res_naive.second.second << endl;

  }

  print_to_plot(rec_runtime, "data_matrix/rec_runtime.dat");
  print_to_plot(row_runtime, "data_matrix/row2_runtime.dat");

  print_to_plot(rec_l2, "data_matrix/rec_l2.dat");
  print_to_plot(row_l2, "data_matrix/row2_l2.dat");

  print_to_plot(rec_l3, "data_matrix/rec_l3.dat");
  print_to_plot(row_l3, "data_matrix/row2_l3.dat");

}

void test_runtime_instructions_row_idx_opt() {

  vector<pair<int, double > > row_std_runtime, row_idx_runtime;
  vector<pair<int, double > > row_std_ins, row_idx_ins;

  int size;

  int testruns = 10;
  int numruns = 50;
  int events[1] = {PAPI_TOT_INS};
  int event_size = 1;

  for (int run=1; run<numruns+1; run++) {

    size = run*15;
    
    cout << "Running test " << run << " of " << numruns << endl;
    
    int m = size;
    int n = size;
    int p = size;
  
    matrix A, B;
  
    resize_matrix(A, m, n);
    resize_matrix(B, n, p);

    pair<double,pair<double,double> > res_row_std = test_row_mult_exclude_preprocess(A, B, events, event_size, testruns);
    row_std_runtime.push_back(make_pair(size, res_row_std.first/1000));
    row_std_ins.push_back(make_pair(size, res_row_std.second.first));

    pair<double,pair<double,double> > res_row_idx = test_row_mult_idx_opt_exclude_preprocess(A, B, events, event_size, testruns);
    row_idx_runtime.push_back(make_pair(size, res_row_idx.first/1000));
    row_idx_ins.push_back(make_pair(size, res_row_idx.second.first));
   
    //cout << res_naive.first/1000 << " " << (long long) res_naive.second.first << " " << (long long) res_naive.second.second << endl;

  }

  print_to_plot(row_std_runtime, "data_matrix/row_std_runtime.dat");
  print_to_plot(row_idx_runtime, "data_matrix/row_idx_runtime.dat");

  print_to_plot(row_std_ins, "data_matrix/row_std_ins.dat");
  print_to_plot(row_idx_ins, "data_matrix/row_idx_ins.dat");

}

void test_runtime_instructions_parallel_row() {

  vector<pair<int, double > > row_std_runtime, row_par_runtime;
  vector<pair<int, double > > row_std_ins, row_par_ins;

  int size;

  int testruns = 10;
  int numruns = 50;
  int events[1] = {PAPI_TOT_INS};
  int event_size = 1;

  for (int run=1; run<numruns+1; run++) {

    size = run*15;
    
    cout << "Running test " << run << " of " << numruns << endl;
    
    int m = size;
    int n = size;
    int p = size;
  
    matrix A, B;
  
    resize_matrix(A, m, n);
    resize_matrix(B, n, p);

    pair<double,pair<double,double> > res_row_std = test_row_mult_idx_opt_exclude_preprocess(A, B, events, event_size, testruns);
    row_std_runtime.push_back(make_pair(size, res_row_std.first/1000));
    row_std_ins.push_back(make_pair(size, res_row_std.second.first));

    pair<double,pair<double,double> > res_row_par = test_row_mult_parallel_exclude_preprocess(A, B, events, event_size, testruns);
    row_par_runtime.push_back(make_pair(size, res_row_par.first/1000));
    row_par_ins.push_back(make_pair(size, res_row_par.second.first));
   
    //cout << res_naive.first/1000 << " " << (long long) res_naive.second.first << " " << (long long) res_naive.second.second << endl;

  }

  print_to_plot(row_std_runtime, "data_matrix/row_std_runtime2.dat");
  print_to_plot(row_par_runtime, "data_matrix/row_par_runtime.dat");

  print_to_plot(row_std_ins, "data_matrix/row_std_ins2.dat");
  print_to_plot(row_par_ins, "data_matrix/row_par_ins.dat");

}

void test_runtime_instructions_parallel_naive() {

  vector<pair<int, double > > naive_std_runtime, naive_par_runtime;
  vector<pair<int, double > > naive_std_ins, naive_par_ins;

  int size;

  int testruns = 10;
  int numruns = 50;
  int events[1] = {PAPI_TOT_INS};
  int event_size = 1;

  for (int run=1; run<numruns+1; run++) {

    size = run*15;
    
    cout << "Running test " << run << " of " << numruns << endl;
    
    int m = size;
    int n = size;
    int p = size;
  
    matrix A, B;
  
    resize_matrix(A, m, n);
    resize_matrix(B, n, p);

    pair<double,pair<double,double> > res_naive_std = test_naive_mult_exclude_preprocess(A, B, events, event_size, testruns);
    naive_std_runtime.push_back(make_pair(size, res_naive_std.first/1000));
    naive_std_ins.push_back(make_pair(size, res_naive_std.second.first));

    pair<double,pair<double,double> > res_naive_par = test_naive_mult_parallel_exclude_preprocess(A, B, events, event_size, testruns);
    naive_par_runtime.push_back(make_pair(size, res_naive_par.first/1000));
    naive_par_ins.push_back(make_pair(size, res_naive_par.second.first));
   
    //cout << res_naive.first/1000 << " " << (long long) res_naive.second.first << " " << (long long) res_naive.second.second << endl;

  }

  print_to_plot(naive_std_runtime, "data_matrix/naive_std_runtime.dat");
  print_to_plot(naive_par_runtime, "data_matrix/naive_par_runtime.dat");

  print_to_plot(naive_std_ins, "data_matrix/naive_std_ins.dat");
  print_to_plot(naive_par_ins, "data_matrix/naive_par_ins.dat");

}

void run_all(int from, int to, const char* dname) {

  int events[2] = {PAPI_L2_TCA,PAPI_L3_TCA};
  test_runtime_mirror(from, to, 10, events, 2, dname, "l2_l3", "l2_tca", "l3_tca");

  int events1[2] = {PAPI_TLB_DM, PAPI_TOT_CYC};
  test_runtime_mirror(from, to, 10, events1, 2, dname, "tlb_cyc", "tlb", "cyc");

  int events3[2] = {PAPI_L2_TCM, PAPI_L2_TCA};
  test_runtime_ratio(from, to, 10, events3, 2, dname, "l2_miss", "ratio");

  int events4[2] = {PAPI_L3_TCM, PAPI_L3_TCA};
  test_runtime_ratio(from, to, 10, events4, 2, dname, "l3_miss", "ratio");
  
  int events5[2] = {PAPI_L2_TCM, PAPI_L3_TCM};
  test_runtime_mirror(from, to, 10, events5, 2, dname, "l2_l3_tcm", "l2_tcm", "l3_tcm");
  
}



int main() {

  //run_all(1, 20, "data_matrix/Test0-20/");
  //run_all(20, 70, "data_matrix/Test20-70/");
  //run_all(100, 200, "data_matrix/Test100-200/");
  //run_all(201, 300, "data_matrix/Test200-300/");
  //run_all(400, 450, "data_matrix/Test400-450/");
  //run_all(501, 600, "data_matrix/Test501-600/");
  //run_all(900, 1100, "data_matrix/Test900-1100/");

  /*
  int from = 100;
  int to = 130;
  int step = 1;
  
  int events[2] = {PAPI_L2_TCA,PAPI_L3_TCA};
  test_runtime_square_mirror(from, to, step, 2, events, 2, "data_matrix/squared4/", "l2_l3", "l2_tca", "l3_tca");

  int events1[2] = {PAPI_TLB_DM, PAPI_TOT_CYC};
  test_runtime_square_mirror(from, to, step, 2, events1, 2, "data_matrix/squared4/", "tlb_cyc", "tlb", "cyc");

  int events3[2] = {PAPI_L2_TCM, PAPI_L2_TCA};
  test_runtime_ratio_square(from, to, step, 2, events3, 2, "data_matrix/squared4/", "l2_miss", "ratio");

  int events4[2] = {PAPI_L3_TCM, PAPI_L3_TCA};
  test_runtime_ratio_square(from, to, step, 2, events4, 2, "data_matrix/squared4/", "l3_miss", "ratio");

  int events5[2] = {PAPI_L2_TCM, PAPI_L3_TCM};
  test_runtime_square_mirror(from, to, step, 2, events5, 2, "data_matrix/squared4/", "l2_l3_tcm", "l2_tcm", "l3_tcm");
  */
  
  /*
  int size = 750;
  
  int testruns = 10;
  int numruns = 50;
  int events[2] = {PAPI_L2_TCA, PAPI_L3_TCA};
  int event_size = 2;

  pair<double,pair<double,double> > res_naive2_par4 = test_naive2_parallel_mult_exclude_preprocess(size, 4, events, event_size, testruns);
  cout << res_naive2_par4.first;
  */
  //int size = 1000;

  /*
  for (int i=1; i<10; i++) {
    int size = i*200;
    cout << "N = " << size << " ----------------------------" << endl;
    cout << "Single:\t\t" << test_runtime_mult_naive2(i*200) / 1000 << endl;
    omp_set_num_threads(2);
    cout << "Par (2 cores):\t" << test_runtime_mult_parallel_naive2(i*200) / 1000 << endl;
    omp_set_num_threads(4);
    cout << "Par (4 cores):\t" << test_runtime_mult_parallel_naive2(i*200) / 1000 << endl;
    cout << endl;
  }
  */
  //for (int i=0; i<3; i++)
  //  cout << test_runtime_mult_parallel_naive2() / 1000 << endl;
  

  
  //test_runtime_l2_l3_accesses();


  
  //int events[2] = {PAPI_L2_TCA,PAPI_L3_TCA};
  //test_runtime_mirror(2, 3, 1, events, 2, "data_matrix/Test/", "l2_l3", "l2_tca", "l3_tca");

  //int events1[2] = {PAPI_TLB_DM, PAPI_TOT_CYC};
  //test_runtime_mirror(2, 3, 1, events1, 2, "data_matrix/Test/", "tlb_cyc", "tlb", "cyc");

  //int events3[2] = {PAPI_L2_TCM, PAPI_L2_TCA};
  //test_runtime_ratio(2, 3, 1, events3, 2, "data_matrix/Test/", "l2_miss", "ratio");

  //int events4[2] = {PAPI_L3_TCM, PAPI_L3_TCA};
  //test_runtime_ratio(2, 3, 1, events4, 2, "data_matrix/Test/", "l3_miss", "ratio");

  
  //return 0;

  //test_runtime_l2_l3_accesses_row_rec();
  //test_runtime_instructions_row_idx_opt();
  //test_runtime_instructions_parallel_row();
  //test_runtime_instructions_parallel_naive();

  /*
  srand (time(NULL));

  for (int i=0; i<matrix_size(A).first; i++)
    for (int j=0; j<matrix_size(A).second; j++) {
      A[i][j] = (i*matrix_size(A).second)+j+1;
    }
  
  for (int i=0; i<matrix_size(B).first; i++)
    for (int j=0; j<matrix_size(B).second; j++) {
      B[i][j] = (i*matrix_size(B).second)+j+1;
    }
  
  //pair<vi,pair<vi,vi> > prv = preprocess_mult_naive(A, B);
  //vector<int> res_row = mult_naive(A, B, prv);

  //print_array_as_matrix(res_row, 4, 4);
  
  matrix res = mult_rec_cutoff(A, B, C, 1, matrix_size(A).first, 1, matrix_size(B).first, 1, matrix_size(B).second, 1);

  // print_matrix(res);
  */

  vector<pair<int, double > > rec_runtime, row_runtime;
  vector<pair<int, double > > rec_l2, row_l2;
  vector<pair<int, double > > rec_l3, row_l3;

  int size = 1 << 11;

  int testruns = 1;
  int numruns = 260;
  int events[2] = {PAPI_TLB_DM, PAPI_L2_TCM};
  int event_size = 2;

  for (int run=20; run<numruns+1; run++) {

    //cout << "Running test " << run << " of " << numruns << endl;
    
    int m = size;
    int n = size;
    int p = size;
  
    matrix A, B;
  
    resize_matrix(A, m, n);
    resize_matrix(B, n, p);

    pair<double,pair<double,double> > res_rec = test_mult_rec_cutoff(A, B, events, event_size, testruns, run);
    rec_runtime.push_back(make_pair(run, res_rec.first/1000));
    rec_l2.push_back(make_pair(run, res_rec.second.first));
    rec_l3.push_back(make_pair(run, res_rec.second.second));
   
    cout << run << ": "  << res_rec.first / 1000 << endl;

  }

  print_to_plot(rec_runtime, "data_matrix/Rec3/rec_runtime_cutoff.dat");
  print_to_plot(rec_l2, "data_matrix/Rec3/rec_tbl_cutoff.dat");
  print_to_plot(rec_l3, "data_matrix/Rec3/rec_l2_tcm_cutoff.dat");

};
