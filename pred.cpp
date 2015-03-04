//TODO: shift left instead of multiply 2
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <algorithm>
#include <vector>
#include <queue>
#include <papi.h>
#include <math.h>
#include <iostream>

#define SIZE 5000000
#define NUM_QUERIES 2
#define MAX_NUM 100000000

using namespace std;

struct Node {
  int data;
  struct Node* left;
  struct Node* right;
};

struct Node* newNode(int data) {
  struct Node* node = (struct Node*) malloc(sizeof(struct Node));
  node->data = data;
  node->left = NULL;
  node->right = NULL;
  return node;
}

struct Node* sorted_array_to_BST(vector<int> &a, int start, int end) {
  if (start > end) return NULL;
  int mid = start + ((end - start) / 2);
  struct Node *root = newNode(a[mid]);
  root->left = sorted_array_to_BST(a, start, mid-1);
  root->right = sorted_array_to_BST(a, mid+1, end);
  return root;
}

//predecessor query using binary search in a sorted array
int pred_sorted_array(int x, vector<int> &sorted_arr) {

  int mi = 0, ma = SIZE-1;
  int mid;
  while (ma >= mi) {
    mid = mi+((ma-mi)/2);
    if (sorted_arr[mid] == x)
      return sorted_arr[mid];
    else if (sorted_arr[mid] < x)
      mi = mid+1;
    else
      ma = mid-1;
  }
  if (sorted_arr[mid] > x && mid == 0) return -1;
  if (sorted_arr[mid] > x) return sorted_arr[mid-1];
  return sorted_arr[mid];
}

//predecessor query using std::lower_bound method in a sorted sorted_array
int pred_lower_bound(int x, vector<int> &sorted_arr) {
  vector<int>::iterator it = lower_bound(sorted_arr.begin(), sorted_arr.end(), x);
  if (*it == x) return x;
  if (it == sorted_arr.begin()) return -1;
  return *(it-1);
}

int tree_minimum(int(*left)(const int, const int),
		 int node, int height) {
  while (--height) {
    node = left(node, height);
  }
  return node;
}

int tree_maximum(int(*right)(const int, const int),
		 int node, int height) {
  while (--height) {
    node = right(node, height);
  }
  return node;
}

int tree_predecessor(int(*left)(const int, const int),
		     int(*right)(const int, const int),
		     vector<int> &arr,
		     int elem,
		     int height,
		     int root) {
  int result = -1;
  int size = (int)arr.size();
  while (height && root < size && root >= 0 && arr[root] != -1) {
    if (arr[root] > elem)
      root = left(root,height);
    else {
      result = arr[root];
      root = right(root, height);
    }
    --height;
  }
  
  return result;
}

inline int bfs_right(int x, int height) {
  return (x+1)<<1;
}

inline int bfs_left(int x, int height) {
  return (x<<1)+1;
}

inline int dfs_left(int x, int height) {
  return x+1;
}

inline int dfs_right(int x, int height) {
  // return x + pow(2, height - 1);
  return x + (1<<(height-1));
}

inline int inorder_left(int x, int height) {
  // return x - pow(2, height - 2);
  return x - (1<<(height-2));
}

inline int inorder_right(int x, int height) {
  // return x + pow(2, height - 2);
  return x + (1<<(height-2));
}

void print_array(vector<int> a) {
  printf("size %d \n", (int)a.size());
  for (size_t i = 0; i < a.size(); i++)
    printf("%d ", a[i]);
  puts("");
}

void preprocess_sorted_array(vector<int> &a) {
  sort(a.begin(), a.end());
}

void preprocess_dfs_array_rec(struct Node *n, int root, int height, vector<int> &res) {
  if (n != NULL)
    res[root] = n->data;
  else return;
  if (dfs_left(root, height) < (int)res.size())
    preprocess_dfs_array_rec(n->left, dfs_left(root, height), height-1, res);
  if (dfs_right(root, height) < (int)res.size())
    preprocess_dfs_array_rec(n->right, dfs_right(root, height), height-1, res);
}

void preprocess_bfs_array_rec(struct Node *n, int root, int height, vector<int> &res) {
  if (n != NULL)
    res[root] = n->data;
  else return;
  if (bfs_left(root, height) < (int)res.size())
    preprocess_bfs_array_rec(n->left, bfs_left(root, height), height-1, res);
  if (bfs_right(root, height) < (int)res.size())
    preprocess_bfs_array_rec(n->right, bfs_right(root, height), height-1, res);
}

void preprocess_inorder_array_rec(struct Node *n, int root, int height, vector<int> &res) {
  if (n != NULL)
    res[root] = n->data;
  else return;
  if (inorder_left(root, height < (int)res.size() && inorder_left(root,height) >= 0))
    preprocess_inorder_array_rec(n->left, inorder_left(root, height), height-1, res);
  if (inorder_right(root, height < (int)res.size() && inorder_right(root,height) >= 0))
    preprocess_inorder_array_rec(n->right, inorder_right(root, height), height-1, res);
}

vector<int> sorted_array_to_dfs_tree(vector<int> &sorted_arr) {
  struct Node *BST = sorted_array_to_BST(sorted_arr, 0, SIZE-1);
  vector<int> result;
  int height = ceil(log2((int)SIZE));
  int size = pow(2,height)-1;
  result.resize(size);
  for (int i = 0; i < size; i++) result[i] = -1;
  preprocess_dfs_array_rec(BST, 0, height, result);
  return result;
}

vector<int> sorted_array_to_bfs_tree(vector<int> &sorted_arr) {
  struct Node *BST = sorted_array_to_BST(sorted_arr, 0, SIZE-1);
  vector<int> result;
  int height = ceil(log2((int)SIZE));
  int size = pow(2,height)-1;
  result.resize(size);
  for (int i = 0; i < size; i++) result[i] = -1;
  preprocess_bfs_array_rec(BST, 0, height, result);
  return result;
}

vector<int> sorted_array_to_inorder_tree(vector<int> &sorted_arr) {
  struct Node *BST = sorted_array_to_BST(sorted_arr, 0, SIZE-1);
  vector<int> result;
  int height = ceil(log2((int)SIZE));
  int size = pow(2,height)-1;
  result.resize(size);
  for (int i = 0; i < size; i++) result[i] = -1;
  preprocess_inorder_array_rec(BST, size/2-1, height, result);
  return result;
}

void test() {

  vector<int> sorted_arr;
  sorted_arr.resize(SIZE);
  //initialize seed
  srand (time(NULL));
  //fill the array with random numbers
  
  for (int i = 0; i < SIZE; i++) {
    int r = rand() % MAX_NUM;
    sorted_arr[i] = r;
  }

  preprocess_sorted_array(sorted_arr);
  // print_array(sorted_arr);
  vector<int> dfs = sorted_array_to_dfs_tree(sorted_arr);
  vector<int> bfs = sorted_array_to_bfs_tree(sorted_arr);
  vector<int> inorder = sorted_array_to_inorder_tree(sorted_arr);
  // print_array(dfs);
  // print_array(bfs);
  // print_array(inorder);
  //vector<int> inorder = sorted_arr;
  int x = rand() % MAX_NUM;
  cout << x << endl;
  int height = ceil(log2((int)SIZE));
  
  int std_pred = pred_lower_bound(x, sorted_arr);
  int dfs_pred = tree_predecessor(dfs_left, dfs_right, dfs, x, height, 0);
  int bfs_pred = tree_predecessor(bfs_left, bfs_right, bfs, x, height, 0);
  int inorder_pred = tree_predecessor(inorder_left, inorder_right, inorder, x, height, inorder.size()/2-1);
  int bin_search = pred_sorted_array(x, sorted_arr);

  printf("%d %d %d %d %d\n", std_pred, dfs_pred, bfs_pred, inorder_pred, bin_search);

}

void tester(int(*left)(const int, const int),
            int(*right)(const int, const int),
            vector<int> &arr,
            int x,
            int height,
            int root,
            vector<int> events) {
  // int numEvents = (int)events.size();
  
  int eventset = PAPI_NULL;
  long long values[1] = {0};

  PAPI_library_init(PAPI_VER_CURRENT);

  PAPI_create_eventset(&eventset);
  PAPI_add_event(eventset, PAPI_L1_DCA);

  PAPI_start(eventset);

  //run predecessor query
  int p = tree_predecessor(left, right, arr, x, height, root);
  cout << p << endl;
  //read & stop
  PAPI_read(eventset, values);
  PAPI_stop(eventset, values);
  //print results
  cout << values[0] << endl;
}

int main() {

  vector<int> events;
  events.push_back(PAPI_L1_DCA);

  // test(); 
  vector<int> sorted_arr;
  sorted_arr.resize(SIZE);
  //initialize seed
  srand (time(NULL));
  //fill the array with random numbers
  
  for (int i = 0; i < SIZE; i++) {
    int r = rand() % MAX_NUM;
    sorted_arr[i] = r;
  }

  preprocess_sorted_array(sorted_arr);
  vector<int> dfs = sorted_array_to_dfs_tree(sorted_arr);
  vector<int> bfs = sorted_array_to_bfs_tree(sorted_arr);
  vector<int> inorder = sorted_array_to_inorder_tree(sorted_arr);
  int x = rand() % MAX_NUM;
  int height = ceil(log2((int)SIZE));

  tester(dfs_left, dfs_right, dfs, x, height, 0, events);
  
  
  return 0;
}

int main2 () {

  int eventset = PAPI_NULL;
  long long values[1] = {0};

  PAPI_library_init(PAPI_VER_CURRENT);

  PAPI_create_eventset(&eventset);
  PAPI_add_event(eventset, PAPI_L1_DCA);

  PAPI_start(eventset);

  vector<int> arr;
  arr.resize(SIZE);
  //initialize seed
  srand (time(NULL));
  //  fill the array with random numbers
  
  for (int i = 0; i < SIZE; i++) {
    int r = rand() % MAX_NUM;
    arr[i] = r;
  }

  PAPI_read(eventset, values);
  printf("number of total instructions:  %lld\n", values[0]);
  PAPI_stop(eventset, values);

  return 0;
}
