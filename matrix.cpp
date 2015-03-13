#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <stdio.h>
#include <algorithm>
#include <papi.h>
#include <fstream>

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

void mult_row_layout(matrix &A, matrix &B, vector<int> &res, double &time, int events[], int event_size, long long &count0, long long &count1) {
  
  clock_t start,end;
  
  vi a, b;
  
  preprocess_matrix_row_layout(A,a);
  preprocess_matrix_col_layout(B,b);

  int a_num_rows = matrix_size(A).first;
  int a_num_cols = matrix_size(A).second;
  int b_num_rows = matrix_size(B).first;
  int b_num_cols = matrix_size(B).second;
  
  res.resize(a_num_rows*b_num_cols,0);

  PAPI_library_init(PAPI_VER_CURRENT);
  long long values[event_size];

  PAPI_start_counters(events, event_size);
  
  start = clock();
  
  for (int k=0; k<b_num_cols;k++) {
    for (int j=0; j<a_num_rows;j++) {
      for (int i=0; i<a_num_cols;i++) {
	res[(k*a_num_rows)+j] = res[(k*a_num_rows)+j] + (a[(j*a_num_cols)+i] * b[(k*b_num_rows)+i]);
      }
    }
  }
  
  PAPI_stop_counters(values, event_size);

  count0 = values[0];
  count1 = values[1];
  
  end = clock();
  time = (end - start) / (double)(CLOCKS_PER_SEC / 1000);
  
}

void mult_col_layout(matrix &A, matrix &B, vector<int> &res, double &time, int events[], int event_size, long long &count0, long long &count1) {
  
  clock_t start,end;
  
  vi a, b;
  
  preprocess_matrix_col_layout(A,a);
  preprocess_matrix_row_layout(B,b);
  
  int a_num_rows = matrix_size(A).first;
  int a_num_cols = matrix_size(A).second;
  int b_num_cols = matrix_size(B).second;
  
  res.resize(a_num_rows*b_num_cols,0);

  PAPI_library_init(PAPI_VER_CURRENT);
  long long values[event_size];

  PAPI_start_counters(events, event_size);

  start = clock();
  
  for (int k=0; k<b_num_cols;k++) {
    for (int j=0; j<a_num_rows;j++) {
      for (int i=0; i<a_num_cols;i++) {
	res[(k*a_num_rows)+j] = res[(k*a_num_rows)+j] + (a[(i*a_num_rows)+j] * b[(i*b_num_cols)+k]);
      }
    }
  }

  PAPI_stop_counters(values, event_size);

  count0 = values[0];
  count1 = values[1];
   
  end = clock();
  time = (end - start) / (double)(CLOCKS_PER_SEC / 1000);
  
}

void mult_naive_layout(matrix &A, matrix &B, vector<int> &res, double &time, int events[], int event_size, long long &count0, long long &count1) {

  clock_t start,end;
  
  vi a, b;
  
  preprocess_matrix_row_layout(A,a);
  preprocess_matrix_row_layout(B,b); 
  
  int a_num_rows = matrix_size(A).first;
  int a_num_cols = matrix_size(A).second;
  int b_num_cols = matrix_size(B).second;
  
  res.resize(a_num_rows*b_num_cols,0);

  PAPI_library_init(PAPI_VER_CURRENT);
  long long values[event_size];

  PAPI_start_counters(events, event_size);
  
  start = clock();

  for (int k=0; k<b_num_cols;k++) {
    for (int j=0; j<a_num_rows;j++) {
      for (int i=0; i<a_num_cols;i++) {
	res[(k*a_num_rows)+j] = res[(k*a_num_rows)+j] + (a[(j*a_num_cols)+i] * b[(i*b_num_cols)+k]);
      }
    }
  }

  end = clock();

  PAPI_stop_counters(values, event_size);

  count0 = values[0];
  count1 = values[1];

  time = (end - start) / (double)(CLOCKS_PER_SEC / 1000);

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

#define MAX_NUM 10000000

void print_to_plot(vector<pair<int,double> > &data, char* fname) {
  ofstream myfile;
  myfile.open(fname);
  cout << "#x\ty" << endl;
  myfile << "#x\ty" << endl;
  for (size_t i = 0; i < data.size(); i++) {
    myfile << data[i].first << "\t" << data[i].second << endl;
    cout << data[i].first << "\t" << data[i].second << endl;
  }
  myfile.close();
}

vector<pair<int, double > > naive_res, col_res, row_res;

void test_running_time(int size) {

  matrix A, B;
  
  resize_matrix(A, size, size);
  resize_matrix(B, size, size);

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

  vector<int> res_naive;
  vector<int> res_col;
  vector<int> res_row;

  double row_time;
  double col_time;
  double naive_time;

  long long naive_count0 = 0;
  long long naive_count1 = 0;

  long long col_count0 = 0;
  long long col_count1 = 0;

  long long row_count0 = 0;
  long long row_count1 = 0;
  
  int events[2] = {PAPI_L2_TCA, PAPI_L2_TCM};
  int event_size = 2;
  
  cout << size << "x" << size << " ----------------------------------------------------" << endl;

  mult_row_layout(A, B, res_row, row_time, events, event_size, row_count0, row_count1);
  cout << "row_mult\t" << "time: " << row_time/1000  << "\t";
  cout << "tca: " << row_count0 << "\ttcm: " << row_count1 << "\t";
  cout << "miss ratio: " << (double) row_count1 / (double) row_count0 << "\t";
  cout << "Cache access / sec: " << (double) row_count0 / row_time << "\t";
  cout << "Cache misses / sec: " << (double) row_count1 / row_time << "\t";
  cout << "sum: " << array_sum(res_row) << "\n";

  mult_naive_layout(A, B, res_naive, naive_time, events, event_size, naive_count0, naive_count1); 
  cout << "naive_mult\t" << "time: " << naive_time/1000  << "\t";
  cout << "tca: " << naive_count0 << "\ttcm: " << naive_count1 << "\t";
  cout << "miss ratio: " << (double) naive_count1 / (double) naive_count0 << "\t";
  cout << "Cache access / sec: " << (double) naive_count0 /  naive_time << "\t";
  cout << "Cache misses / sec: " << (double) naive_count1 /  naive_time << "\t";
  cout << "sum: " << array_sum(res_naive) << endl;

  mult_col_layout(A, B, res_col, col_time, events, event_size, col_count0, col_count1);
  cout << "col_mult\t" << "time: " << col_time/1000  << "\t";
  cout << "tca: " << col_count0 << "\ttcm: " << col_count1 << "\t";
  cout << "miss ratio: " << (double) col_count1 / (double) col_count0 << "\t";
  cout << "Cache access / sec: " << (double) col_count0 / col_time << "\t";
  cout << "Cache misses / sec: " << (double) col_count1 / col_time << "\t";
  cout << "sum: " << array_sum(res_col) << endl;

  naive_res.push_back(make_pair(size, (double) row_count1 / row_time));
  col_res.push_back(make_pair(size, (double) col_count1 / row_time));
  row_res.push_back(make_pair(size, (double) row_count1 / row_time));
  
  print_to_plot(naive_res, "Naive.dat");
  print_to_plot(col_res, "Col.dat");
  print_to_plot(row_res, "Row.dat");
  
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
  
}

void test_mult_rec() {

   matrix A, B, C;
  
  resize_matrix(A,4,5);
  resize_matrix(B,5,4);
  resize_matrix(C,4,4);

  for (int i=0; i<matrix_size(A).first; i++)
    for (int j=0; j<matrix_size(A).second; j++)
      A[i][j] = (i*matrix_size(A).second)+j+1;

  for (int i=0; i<matrix_size(B).first; i++)
    for (int j=0; j<matrix_size(B).second; j++)
      B[i][j] = (i*matrix_size(B).second)+j+1;

  matrix res = mult_rec(A, B, C, 1, matrix_size(A).first, 1, matrix_size(B).first, 1, matrix_size(B).second);
  
  print_matrix(res);

}

int main() {

  //test_mult_rec();
  test_running_time(200);

  return 0;
  
};



