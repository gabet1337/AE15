#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <algorithm>
#include <vector>
#include <queue>
#include <papi.h>
#include <math.h>
#include <iostream>
#include <fstream>

//2e12*4
//#define SIZE 16383
//2e16*4
//#define SIZE 262143
//2e19*4

#define NUM_QUERIES 1e6
// #define NUM_QUERIES 1
#define MAX_NUM 10000000
#define NUM_EXPERIMENTS 42

using namespace std;

void clear_cache() {
  ofstream ofs("/proc/sys/vm/drop_caches");
  ofs << "3" << endl;
}

void print_to_plot(vector<pair<int,long long> > &data) {
  cout << "#x\ty" << endl;
  for (size_t i = 0; i < data.size(); i++)
    cout << data[i].first << "\t" << data[i].second << endl;
}

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

  int mi = 0, ma = (int)sorted_arr.size()-1;
  int mid = 0;
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

int bfs_right(int x, int height) {
  return (x+1)<<1;
}

int bfs_left(int x, int height) {
  return (x<<1)+1;
}

int dfs_left(int x, int height) {
  return x+1;
}

int dfs_right(int x, int height) {
  // return x + pow(2, height - 1);
  return x + (1<<(height-1));
}

int inorder_left(int x, int height) {
  // return x - pow(2, height - 2);
  return x - (1<<(height-2));
}

int inorder_right(int x, int height) {
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

vector<int> sorted_array_to_dfs_tree(vector<int> &sorted_arr, size_t s) {
  struct Node *BST = sorted_array_to_BST(sorted_arr, 0, s-1);
  vector<int> result;
  int height = ceil(log2((int)s+1));
  int size = pow(2,height+1)-1;
  result.resize(size);
  for (int i = 0; i < size; i++) result[i] = -1;
  preprocess_dfs_array_rec(BST, 0, height+1, result);
  return result;
}

vector<int> sorted_array_to_bfs_tree(vector<int> &sorted_arr, size_t s) {
  struct Node *BST = sorted_array_to_BST(sorted_arr, 0, s-1);
  vector<int> result;
  int height = ceil(log2((int)s+1));
  int size = pow(2,height+1)-1;
  result.resize(size);
  for (int i = 0; i < size; i++) result[i] = -1;
  preprocess_bfs_array_rec(BST, 0, height+1, result);
  return result;
}

vector<int> sorted_array_to_inorder_tree(vector<int> &sorted_arr, size_t s) {
  struct Node *BST = sorted_array_to_BST(sorted_arr, 0, s-1);
  vector<int> result;
  int height = ceil(log2((int)s+1));
  int size = pow(2,height+1)-1;
  result.resize(size);
  for (int i = 0; i < size; i++) result[i] = -1;
  preprocess_inorder_array_rec(BST, size/2-1, height+1, result);
  return result;
}

void test2() {
  for (int expo = 3; expo <= 22; expo++) {
    vector<int> sorted_arr;
    int s = pow(2,expo)-1;

    sorted_arr.resize(s);
    //initialize seed
    srand (time(NULL));
    //fill the array with random numbers
  
    for (int i = 0; i < s; i++) {
      int r = rand() % MAX_NUM;
      // cout << r << endl;
      sorted_arr[i] = r;
    }

    preprocess_sorted_array(sorted_arr);
    // print_array(sorted_arr);
    //vector<int> dfs = sorted_array_to_dfs_tree(sorted_arr,s);
    //vector<int> bfs = sorted_array_to_bfs_tree(sorted_arr,s);
    vector<int> inorder = sorted_array_to_inorder_tree(sorted_arr,s);
    // print_array(inorder);
    // print_array(dfs);
    // print_array(bfs);
    // print_array(inorder);
    //vector<int> inorder = sorted_arr;
    PAPI_library_init(PAPI_VER_CURRENT);
    int eventset = PAPI_NULL;
    PAPI_create_eventset(&eventset);
    PAPI_add_event(eventset, PAPI_L2_TCM);
    long long values[1] = {0};
    PAPI_start(eventset);
    for (int i = 0; i < 1e6; i++) {
      int x = rand() % MAX_NUM;
      int height = ceil(log2((int)s));
  
      //int std_pred = pred_lower_bound(x, sorted_arr);
      //int dfs_pred = tree_predecessor(dfs_left, dfs_right, dfs, x, height, 0);
      //int bfs_pred = tree_predecessor(bfs_left, bfs_right, bfs, x, height, 0);
      int inorder_pred = tree_predecessor(inorder_left, inorder_right, inorder, x, height, inorder.size()/2-1);
      inorder_pred++;
      // printf("predecessor of %d is %d\n",x, inorder_pred); 
      //int bin_search = pred_sorted_array(x, sorted_arr);

      //printf("%d %d %d %d %d\n", std_pred, dfs_pred, bfs_pred, inorder_pred, bin_search);
    }
    PAPI_read(eventset,values);
    cout << expo << " " <<  values[0] << endl;
    PAPI_stop(eventset, values);

  }
}

void test() {

  vector<int> sorted_arr;
  int s = (1<<20)-1;
  sorted_arr.resize(s);
  //initialize seed
  srand (time(NULL));
  //fill the array with random numbers
  
  for (int i = 0; i < s; i++) {
    int r = rand() % MAX_NUM;
    sorted_arr[i] = r;
  }

  preprocess_sorted_array(sorted_arr);
  // print_array(sorted_arr);
  vector<int> dfs = sorted_array_to_dfs_tree(sorted_arr,s);
  vector<int> bfs = sorted_array_to_bfs_tree(sorted_arr,s);
  vector<int> inorder = sorted_array_to_inorder_tree(sorted_arr,s);
  // print_array(dfs);
  // print_array(bfs);
  // print_array(inorder);
  //vector<int> inorder = sorted_arr;
  PAPI_library_init(PAPI_VER_CURRENT);
  int eventset = PAPI_NULL;
  PAPI_create_eventset(&eventset);
  PAPI_add_event(eventset, PAPI_L1_TCM);
  long long values[1] = {0};
  PAPI_start(eventset);
  for (int i = 0; i < 100; i++) {
    int x = rand() % MAX_NUM;
    cout << x << endl;
    int height = ceil(log2((int)s));
  
    int std_pred = pred_lower_bound(x, sorted_arr);
    int dfs_pred = tree_predecessor(dfs_left, dfs_right, dfs, x, height, 0);
    int bfs_pred = tree_predecessor(bfs_left, bfs_right, bfs, x, height, 0);
    int inorder_pred = tree_predecessor(inorder_left, inorder_right, inorder, x, height, inorder.size()/2-1);
    int bin_search = pred_sorted_array(x, sorted_arr);

    printf("%d %d %d %d %d\n", std_pred, dfs_pred, bfs_pred, inorder_pred, bin_search);
  }
  PAPI_read(eventset,values);
  cout << values[0] << endl;
  PAPI_stop(eventset, values);
}

void test_bfs() {
  PAPI_library_init(PAPI_VER_CURRENT);
  for (size_t i = 2; i < 25; i++) {
    size_t s = (1<<i)-1;
    vector<int> sorted_arr;
    sorted_arr.resize(s);
    for (size_t k = 0; k < s; k++) {
      sorted_arr[k] = rand() % MAX_NUM;
    }
    preprocess_sorted_array(sorted_arr);
    vector<int> bfs = sorted_array_to_bfs_tree(sorted_arr, s);
    int height = ceil(log2((int)s));
    
    int eventset = PAPI_NULL;
    PAPI_create_eventset(&eventset);
    PAPI_add_event(eventset, PAPI_L2_TCM);
    long long values[1] = {(long long) 0};
    PAPI_start(eventset);
    srand(time(NULL));
    for (size_t j = 0; j < 1e6; j++) {
      int p = tree_predecessor(bfs_left, bfs_right, bfs, rand() % MAX_NUM, height, 0);
      cout << p << endl;
    }
    PAPI_read(eventset, values);
    cout << i << "\t" << values[0] << endl;
    PAPI_stop(eventset, values);

  }
}

long long tester(int(*left)(const int, const int),
                 int(*right)(const int, const int),
                 vector<int> &arr,
                 int height,
                 int root) {

  PAPI_library_init(PAPI_VER_CURRENT);

  int evts[2] = {PAPI_TOT_INS, PAPI_L2_TCA};
  long long values[2];
  
  long long average = 0;

  for (int ex = 0; ex < NUM_EXPERIMENTS; ex++) {
    clear_cache();

    PAPI_start_counters(evts, 2);

    //run predecessor query
    for (size_t i = 0; i < NUM_QUERIES; i++) {
      int r = rand() % MAX_NUM;
      int p = tree_predecessor(left, right, arr, r, height, root);
      p++;
    }
    //read & stop
    PAPI_stop_counters(values, 2);
    // cout << values[0] << endl;
    average += values[1];

  }

  return average/NUM_EXPERIMENTS;
}

int main() {
  srand (time(NULL));
  //test_bfs();
  //test();
  //return 0;

  vector<pair<int,long long> > inorder_res, dfs_res, bfs_res;
  //initialize seed
  for (int ex = 9; ex <= 14; ex++) {
    size_t s = pow(2,ex)-1;
    cout << "Running experiment on input size: " << s << endl;

    vector<int> sorted_arr;
    sorted_arr.resize(s);
    
    //fill the array with random numbers
    for (size_t i = 0; i < s; i++) {
      int r = rand() % MAX_NUM;
      sorted_arr[i] = r;
    }

    preprocess_sorted_array(sorted_arr);

    vector<int> dfs = sorted_array_to_dfs_tree(sorted_arr,s);
    vector<int> bfs = sorted_array_to_bfs_tree(sorted_arr,s);
    vector<int> inorder = sorted_array_to_inorder_tree(sorted_arr,s);
    // print_array(dfs);
    // print_array(bfs);
    // print_array(inorder);
    int inorder_height = ceil(log2((int)inorder.size()));
    int bfs_height = ceil(log2((int)bfs.size()));
    int dfs_height = ceil(log2((int)dfs.size()));
    int inorder_root = inorder.size()/2-1;

    bfs_res.push_back(make_pair(ex, tester(bfs_left, bfs_right, bfs,
                                           bfs_height, 0)));

    
    dfs_res.push_back(make_pair(ex, tester(dfs_left, dfs_right, dfs,
                                           dfs_height,0)));
    
    inorder_res.push_back(make_pair
                          (ex, tester(inorder_left, inorder_right,
                                      inorder, inorder_height, inorder_root)));
  } 
  cout << "INORDER" << endl;
  print_to_plot(inorder_res);
  cout << "BFS" << endl;
  print_to_plot(bfs_res);
  cout << "DFS" << endl;
  print_to_plot(dfs_res);

  return 0;
}
