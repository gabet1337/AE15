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
#define MAX_NUM 100000000
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
  int end = (int)sorted_arr.size()-1;
  struct Node *BST = sorted_array_to_BST(sorted_arr, 0, end);
  vector<int> result;
  int height = ceil(log2((int)end+1));
  int size = pow(2,height)-1;
  result.resize(size);
  for (int i = 0; i < size; i++) result[i] = -1;
  preprocess_dfs_array_rec(BST, 0, height, result);
  return result;
}

vector<int> sorted_array_to_bfs_tree(vector<int> &sorted_arr) {
  int end = (int)sorted_arr.size()-1;
  struct Node *BST = sorted_array_to_BST(sorted_arr, 0, end);
  vector<int> result;
  int height = ceil(log2((int)end+1));
  int size = pow(2,height)-1;
  result.resize(size);
  for (int i = 0; i < size; i++) result[i] = -1;
  preprocess_bfs_array_rec(BST, 0, height, result);
  return result;
}

vector<int> sorted_array_to_inorder_tree(vector<int> &sorted_arr) {
  int end = (int)sorted_arr.size()-1;
  struct Node *BST = sorted_array_to_BST(sorted_arr, 0, end);
  vector<int> result;
  int height = ceil(log2((int)end+1));
  int size = pow(2,height)-1;
  result.resize(size);
  for (int i = 0; i < size; i++) result[i] = -1;
  preprocess_inorder_array_rec(BST, size/2-1, height, result);
  return result;
}

// void test() {

//   vector<int> sorted_arr;
//   sorted_arr.resize(SIZE);
//   //initialize seed
//   srand (time(NULL));
//   //fill the array with random numbers
  
//   for (int i = 0; i < SIZE; i++) {
//     int r = rand() % MAX_NUM;
//     sorted_arr[i] = r;
//   }

//   preprocess_sorted_array(sorted_arr);
//   // print_array(sorted_arr);
//   vector<int> dfs = sorted_array_to_dfs_tree(sorted_arr);
//   vector<int> bfs = sorted_array_to_bfs_tree(sorted_arr);
//   vector<int> inorder = sorted_array_to_inorder_tree(sorted_arr);
//   // print_array(dfs);
//   // print_array(bfs);
//   // print_array(inorder);
//   //vector<int> inorder = sorted_arr;
//   int x = rand() % MAX_NUM;
//   cout << x << endl;
//   int height = ceil(log2((int)SIZE));
  
//   int std_pred = pred_lower_bound(x, sorted_arr);
//   int dfs_pred = tree_predecessor(dfs_left, dfs_right, dfs, x, height, 0);
//   int bfs_pred = tree_predecessor(bfs_left, bfs_right, bfs, x, height, 0);
//   int inorder_pred = tree_predecessor(inorder_left, inorder_right, inorder, x, height, inorder.size()/2-1);
//   int bin_search = pred_sorted_array(x, sorted_arr);

//   printf("%d %d %d %d %d\n", std_pred, dfs_pred, bfs_pred, inorder_pred, bin_search);

// }

long long tester(int(*left)(const int, const int),
            int(*right)(const int, const int),
            vector<int> &arr,
            int height,
            int root,
            vector<int> &events,
            int num_queries) {
  // int numEvents = (int)events.size();
  
  int eventset = PAPI_NULL;
  long long values[1] = {0};

  PAPI_library_init(PAPI_VER_CURRENT);

  PAPI_create_eventset(&eventset);
  for (size_t i = 0; i < events.size(); i++)
    PAPI_add_event(eventset, events[i]);

  long long average = 0;
  for (int ex = 0; ex < NUM_EXPERIMENTS; ex++) {
    clear_cache();
    PAPI_start(eventset);

    //run predecessor query
    for (int i = 0; i < num_queries; i++) {
      int r = rand() % MAX_NUM;
      int p = tree_predecessor(left, right, arr, r, height, root);
      p++;
    }
    //read & stop
    PAPI_read(eventset, values);
    PAPI_stop(eventset, values);
    //print results
    average += values[0];
  }
  return average/NUM_EXPERIMENTS;
}

int main() {

  vector<int> events;
  events.push_back(PAPI_L1_TCM);

  vector<pair<int,long long> > inorder_res, dfs_res, bfs_res;

  for (int ex = 2; ex < 24; ex++) {
    cout << "Running experiment on input size: " << (1<<ex)-1 << endl;
    int size = (1<<ex)-1;
    vector<int> sorted_arr;
    sorted_arr.resize(size);
    //initialize seed
    srand (time(NULL));
    //fill the array with random numbers

    for (int i = 0; i < size; i++) {
      int r = rand() % MAX_NUM;
      sorted_arr[i] = r;
    }

    preprocess_sorted_array(sorted_arr);

    vector<int> inorder = sorted_array_to_inorder_tree(sorted_arr);

    int height = ceil(log2((int)size));
    int inorder_root = inorder.size()/2-1;
    inorder_res.push_back(make_pair(ex, tester(inorder_left, inorder_right,
                                               inorder, height, inorder_root,
                                               events, NUM_QUERIES)));

    vector<int> dfs = sorted_array_to_dfs_tree(sorted_arr);
    dfs_res.push_back(make_pair(ex, tester(dfs_left, dfs_right, dfs,
                                           height, 0, events, NUM_QUERIES)));
  
    vector<int> bfs = sorted_array_to_bfs_tree(sorted_arr);
    bfs_res.push_back(make_pair(ex, tester(bfs_left, bfs_right, bfs,
                                           height, 0, events, NUM_QUERIES)));
  }
  cout << "INORDER" << endl;
  print_to_plot(inorder_res);
  cout << "BFS" << endl;
  print_to_plot(bfs_res);
  cout << "DFS" << endl;
  print_to_plot(dfs_res);
  return 0;
}
