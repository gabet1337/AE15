#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <stdio.h>
#include <algorithm>

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

void mult_row_layout(matrix &A, matrix &B, vector<int> &res, double &time) {

  clock_t start,end;
  
  vi a, b;

  preprocess_matrix_row_layout(A,a);
  preprocess_matrix_col_layout(B,b);

  int a_num_rows = matrix_size(A).first;
  int a_num_cols = matrix_size(A).second;
  int b_num_rows = matrix_size(B).first;
  int b_num_cols = matrix_size(B).second;
  
  res.resize(a_num_rows*b_num_cols,0);

  start = clock();
  
  for (int k=0; k<b_num_cols;k++) {
    for (int j=0; j<a_num_rows;j++) {
      for (int i=0; i<a_num_cols;i++) {
	res[(k*a_num_rows)+j] = res[(k*a_num_rows)+j] + (a[(j*a_num_cols)+i] * b[(k*b_num_rows)+i]);
      }
    }
  }

  end = clock();
  time = (end - start) / (double)(CLOCKS_PER_SEC / 1000);
  
}

void mult_col_layout(matrix &A, matrix &B, vector<int> &res, double &time) {

  clock_t start,end;
  
  vi a, b;
  
  preprocess_matrix_col_layout(A,a);
  preprocess_matrix_row_layout(B,b);

  int a_num_rows = matrix_size(A).first;
  int a_num_cols = matrix_size(A).second;
  //int b_num_rows = matrix_size(B).first;
  int b_num_cols = matrix_size(B).second;
  
  res.resize(a_num_rows*b_num_cols,0);

  start = clock();
  
  for (int k=0; k<b_num_cols;k++) {
    for (int j=0; j<a_num_rows;j++) {
      for (int i=0; i<a_num_cols;i++) {
	res[(k*a_num_rows)+j] = res[(k*a_num_rows)+j] + (a[(i*a_num_rows)+j] * b[(i*b_num_cols)+k]);
      }
    }
  }

  end = clock();
  time = (end - start) / (double)(CLOCKS_PER_SEC / 1000);

}

void mult_naive_layout(matrix &A, matrix &B, vector<int> &res, double &time) {

  clock_t start,end;
  
  vi a, b;

  preprocess_matrix_row_layout(A,a);
  preprocess_matrix_row_layout(B,b); 

  int a_num_rows = matrix_size(A).first;
  int a_num_cols = matrix_size(A).second;
  //int b_num_rows = matrix_size(B).first;
  int b_num_cols = matrix_size(B).second;

  res.resize(a_num_rows*b_num_cols,0);

  start = clock();
  
  for (int k=0; k<b_num_cols;k++) {
    for (int j=0; j<a_num_rows;j++) {
      for (int i=0; i<a_num_cols;i++) {
	res[(k*a_num_rows)+j] = res[(k*a_num_rows)+j] + (a[(j*a_num_cols)+i] * b[(i*b_num_cols)+k]);
      }
    }
  }

  end = clock();
  time = (end - start) / (double)(CLOCKS_PER_SEC / 1000);

}

int array_sum(vector<int> &a) {
  int res = 0;
  for (size_t i=0; i<a.size(); i++)
    res += a[i];

  return res;
}

void mult2_rec(matrix A, matrix B, matrix C) {

  
  
}

void mult_naive_layout(matrix &A, int A_row_idx, int A_rows, int A_col_idx, int A_cols, matrix &B, int B_row_idx, int B_rows, int B_col_idx, int B_cols, matrix &C) {

  //clock_t start,end;
  
  //vi a, b;

  //preprocess_matrix_row_layout(A,a);
  //preprocess_matrix_row_layout(B,b); 

  int a_num_rows = A_rows; // matrix_size(A).first;
  int a_num_cols = A_cols; // matrix_size(A).second;
  //int b_num_rows = matrix_size(B).first;
  int b_num_cols = B_cols; //matrix_size(B).second;

  //res.resize(a_num_rows*b_num_cols,0);

  //start = clock();
  
  for (int k=0; k<b_num_cols;k++) {
    for (int j=0; j<a_num_rows;j++) {
      for (int i=0; i<a_num_cols;i++) {
		C[k+B_col_idx][j+A_row_idx] = C[k+B_col_idx][j+A_row_idx] + (A[j+A_row_idx][i+A_col_idx] * B[i+B_row_idx][k+B_col_idx]);
      }
    }
  }

  //end = clock();
  //time = (end - start) / (double)(CLOCKS_PER_SEC / 1000);

}

void mult_rec(matrix &A, matrix &B, matrix &C, int m_start, int m_end, int n_start, int n_end, int p_start,  int p_end) {
  
  // Invariant: A is mxn, B is nxp, C is mxp
  // M[row][column]

  int m = m_end - m_start + 1;
  int n = n_end - n_start + 1;
  int p = p_end - p_start + 1;
   
  if ((m == 1) && (n == 1) && (p == 1)) {
    C[m_start-1][p_start-1] = C[m_start-1][p_start-1] + (A[m_start-1][n_start-1] * B[n_start-1][p_start-1]);
  } else if (m >= max(n, p)) {
    cout << "a" << endl;
    int m_mid = m/2; // split m
    
	//mult_naive_layout(matrix &A, int A_row_idx, int A_rows, int A_col_idx, int A_cols, matrix &B, int B_row_idx, int B_rows, int B_col_idx, int B_cols, matrix &C)

	mult_rec(A, B, C, m_start, m_end-m_mid, n_start, n_end, p_start, p_end);
	mult_naive_layout(A, m_start, m_end-m_mid+1, n_start, n, B, n_start, n, p_start, p, C);

	mult_rec(A, B, C,m_end-m_mid+1, m_end, n_start, n_end, p_start, p_end);



  } else if (n >= max(m, p)) {
    cout << "b" << endl;
    int n_mid = n/2; // split n
    mult_rec(A, B, C, m_start, m_end, n_start, n_end-n_mid, p_start, p_end);
    mult_rec(A, B, C, m_start, m_end, n_end-n_mid+1, n_end, p_start, p_end);
  } else if (p >= max(m, n)) {
    cout << "b" << endl;
    int p_mid = p/2; //split p
    mult_rec(A, B, C, m_start, m_end, n_start, n_end, p_start, p_end-p_mid);
    mult_rec(A, B, C, m_start, m_end, n_start, n_end, p_end-p_mid+1, p_end);
  }
 
}
      
int main() {

  matrix A, B;

  resize_matrix(A,4,4);
  resize_matrix(B,4,4);

  for (int i=0; i<matrix_size(A).first; i++)
    for (int j=0; j<matrix_size(A).second; j++)
      A[i][j] = (i*matrix_size(A).second)+j+1;

  for (int i=0; i<matrix_size(B).first; i++)
    for (int j=0; j<matrix_size(B).second; j++)
      B[i][j] = (i*matrix_size(B).second)+j+1;

  // vector<int> res_naive;
  // vector<int> res_col;
  // vector<int> res_row;

  // //double row_time;
  // double col_time;
  // double naive_time;
  
  //mult_naive_layout(A, B, res_naive, naive_time); 
  //cout << "naive_mult\t" << "time: " << naive_time/1000  << "\t" << "sum:\t" << array_sum(res_naive) << "\n";

  //mult_col_layout(A, B, res_col, col_time);
  //cout << "col_mult\t" << "time: " << col_time/1000  << "\t" << "sum:\t" << array_sum(res_col) << "\n";

  //mult_row_layout(A, B, res_row, row_time);
  //cout << "row_mult\t" << "time: " << row_time/1000  << "\t" << "sum:\t" << array_sum(res_row) << "\n";

  matrix C;
  resize_matrix(C,4,4);
  
  mult_rec(A, B, C, 1, 4, 1, 4, 4, 4);
  
  //print_array_as_matrix(res_row, matrix_size(A).first, matrix_size(B).second);
  
  return 0;
};

