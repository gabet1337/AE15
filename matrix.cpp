#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <stdio.h>
#include <algorithm>
#include <papi.h>
#include <fstream>
#include <omp.h>

using namespace std;

typedef vector<int> vi;
typedef vector< vi > matrix;
typedef pair<int,int> ii;


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

#define MAX_NUM 10

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

pair<double, pair<double,double> > test_row_mult_exclude_preprocess(matrix &A, matrix &B, int events[], int event_size, int numruns) {
  
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
    vector<int> res_row = mult_row(A, B, prv);
  
    end = clock();

    PAPI_stop_counters(values, event_size);
  
    count0 += values[0];
  
    if (event_size == 2)
      count1 += values[1];
  
    timer += (end - start) / (double)(CLOCKS_PER_SEC / 1000);

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
    vector<int> res = mult_col(A, B, prv);
  
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

pair<double, pair<double,double> > test_naive_mult_exclude_preprocess(matrix &A, matrix &B, int events[], int event_size, int numruns) {
  
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
    vector<int> res = mult_naive(A, B, prv);
  
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

  vector<pair<int, double > > naive_runtime, col_runtime, row_runtime;
  vector<pair<int, double > > naive_l2, col_l2, row_l2;
  vector<pair<int, double > > naive_l3, col_l3, row_l3;

  int size;

  int testruns = 10;
  int numruns = 50;
  int events[2] = {PAPI_L2_TCA, PAPI_L3_TCA};
  int event_size = 2;

  for (int run=1; run<numruns+1; run++) {

    size = run*15;

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

    if (run <= 36) {
      pair<double,pair<double,double> > res_col = test_col_mult_exclude_preprocess(A, B, events, event_size, testruns);
      col_runtime.push_back(make_pair(size, res_col.first/1000));
      col_l2.push_back(make_pair(size, res_col.second.first));
      col_l3.push_back(make_pair(size, res_col.second.second));
    }

    pair<double,pair<double,double> > res_naive = test_naive_mult_exclude_preprocess(A, B, events, event_size, testruns);
    naive_runtime.push_back(make_pair(size, res_naive.first/1000));
    naive_l2.push_back(make_pair(size, res_naive.second.first));
    naive_l3.push_back(make_pair(size, res_naive.second.second));    

    //cout << res_naive.first/1000 << " " << (long long) res_naive.second.first << " " << (long long) res_naive.second.second << endl;

  }

  print_to_plot(naive_runtime, "data_matrix/naive_runtime.dat");
  print_to_plot(col_runtime, "data_matrix/col_runtime.dat");
  print_to_plot(row_runtime, "data_matrix/row_runtime.dat");

  print_to_plot(naive_l2, "data_matrix/naive_l2.dat");
  print_to_plot(col_l2, "data_matrix/col_l2.dat");
  print_to_plot(row_l2, "data_matrix/row_l2.dat");

  print_to_plot(naive_l3, "data_matrix/naive_l3.dat");
  print_to_plot(col_l3, "data_matrix/col_l3.dat");
  print_to_plot(row_l3, "data_matrix/row_l3.dat");

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



int main() {


  clock_t start,end;
  double timer = 0;

  const int size = 4;
  
  int a[size][size];
  int b[size][size];
  int c[size][size];
  
  // Initialize buffers.
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      a[i][j] = i*size+j+1; //  (float)i + j;
      b[i][j] = i*size+j+1; //(float)i - j;
      c[i][j] = 0;
    }
  }
  
  start = clock();

  // Compute matrix multiplication.
  // C <- C + A x B
  //#pragma omp parallel for default(none) shared(a,b,c)
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      for (int k = 0; k < size; ++k) {
	c[i][j] += a[i][k] * b[k][j];
      }
    }
  }
  
  end = clock();

  timer = (end - start) / (double)(CLOCKS_PER_SEC / 1000);

  for (size_t i = 0; i < size; ++i) {
    for (size_t j = 0; j < size; ++j) {
      cout << c[i][j] << '\t';
    }
    cout << endl;
  }
  
  
  cout << timer;
    
  return 0;
    
  //test_runtime_l2_l3_accesses();
  //test_runtime_l2_l3_accesses_row_rec();
  //test_runtime_instructions_row_idx_opt();
  //test_runtime_instructions_parallel_row();
  //test_runtime_instructions_parallel_naive();

  
  /*
  matrix A, B, C;

  int size = 4;
  resize_matrix(A, size, size);
  resize_matrix(B, size, size);
  resize_matrix(C, size, size);
  
  srand (time(NULL));

  for (int i=0; i<matrix_size(A).first; i++)
    for (int j=0; j<matrix_size(A).second; j++) {
      A[i][j] = (i*matrix_size(A).second)+j+1;
    }
  
  for (int i=0; i<matrix_size(B).first; i++)
    for (int j=0; j<matrix_size(B).second; j++) {
      B[i][j] = (i*matrix_size(B).second)+j+1;
    }

  
  pair<vi,pair<vi,vi> > prv = preprocess_mult_naive(A, B);
  vector<int> res_row = mult_naive(A, B, prv);

  print_array_as_matrix(res_row, 4, 4);
  
  //matrix res = mult_rec(A, B, C, 1, matrix_size(A).first, 1, matrix_size(B).first, 1, matrix_size(B).second);

  //print_matrix(res);
  */
  
  return 0;
  
};
