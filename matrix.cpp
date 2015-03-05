#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <stdio.h>

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

void mult_row_layout(int a_num_rows, int a_num_cols, vector<int> &a, int b_num_rows, int b_num_cols, vector<int> &b) {

  vector<int> res;
  res.resize(a_num_rows*b_num_cols,0);

  for (int k=0; k<b_num_cols;k++) {
    for (int j=0; j<a_num_rows;j++) {
      for (int i=0; i<a_num_cols;i++) {
	res[(k*a_num_rows)+j] = res[(k*a_num_rows)+j] + (a[(j*a_num_cols)+i] * b[(k*b_num_rows)+i]);
      }
    }
  }

  cout << "\nResult:\n";
  print_array_as_matrix(res,a_num_rows,b_num_cols);
  
}

void test_row_mult() {

  matrix A, B;

  resize_matrix(A,5,8);
  resize_matrix(B,8,3);

  for (int i=0; i<matrix_size(A).first; i++)
    for (int j=0; j<matrix_size(A).second; j++)
      A[i][j] = (i*matrix_size(A).second)+j+1;

  for (int i=0; i<matrix_size(B).first; i++)
    for (int j=0; j<matrix_size(B).second; j++)
      B[i][j] = (i*matrix_size(B).second)+j+1;
  
  vi a, b;
  
  preprocess_matrix_row_layout(A,a);
  preprocess_matrix_col_layout(B,b);

  cout << "Matrix A:\n";
  print_matrix(A);
 
  cout << "\nMatrix A (row layout)\n";
  print_array(a);

  cout << "\n\nMatrix B:\n";
  print_matrix(B);

  cout << "\nMatrix B (row layout)\n";
  print_array(b);
  
  cout << "\n";

  mult_row_layout(matrix_size(A).first,matrix_size(A).second,a,matrix_size(B).first,matrix_size(B).second,b);
  
}

int main() {

  test_row_mult();
  
  return 0;
};

