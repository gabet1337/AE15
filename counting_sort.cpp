#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

typedef vector<int> vi;

//counting sort as in CLRS page 195
//A is input
//B is output
//k is the range of the integers
void clrs_counting_sort(vi &A, vi &B, int k) {
  vi C(k, 0);
  for (size_t j = 0; j < A.size(); j++)
    C[A[j]]++;
  for (int i = 1; i < k; i++)
    C[i] += C[i-1];
  for (int j = (int)A.size()-1; j >= 0; j--) {
    B[C[A[j]]] = A[j];
    C[A[j]]--;
  }
}

void test_clrs_counting_sort() {
  int num_to_sort = 1000;
  int k = 100;
  vi A,B;

  for (int i = 0; i < num_to_sort; i++)
    A.push_back(rand()%k);

  A.resize(num_to_sort);
  B.resize(num_to_sort);
  clrs_counting_sort(A,B,k);

  if (is_sorted(B.begin(), B.end()))
    cout << "CLRS COUNTING SORT: SUCCESS!" << endl;
}

void article_counting_sort(vi &A, vi &B, int k) {
  size_t n = A.size();
  long long storage[n*k];
  long long next[k];
  for (int i = 0; i < k; i++) next[i] = i*n;
  for (size_t i = 0; i < A.size(); i++) {
    storage[next[A[i]]] = A[i];
    next[A[i]]++;
  }
  int cur = 0;
  for (int i = 0; i < k; i++)
    for (long long j = i*n; j < next[i]; j++)
      B[cur++] = storage[j];
}

void test_article_counting_sort() {
  int num_to_sort = 1000;
  int k = 1000;
  vi A,B;

  for (int i = 0; i < num_to_sort; i++)
    A.push_back(rand()%k);

  A.resize(num_to_sort);
  B.resize(num_to_sort);
  article_counting_sort(A,B,k);

  if (is_sorted(B.begin(), B.end()))
    cout << "ARTICLE COUNTING SORT: SUCCESS!" << endl;

}


int main() {
  srand(time(NULL));
  test_clrs_counting_sort();
  test_article_counting_sort();
  return 0;
}
