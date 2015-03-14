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
#include <string>
#include <cstring>
#include <chrono>
#include <random>

#define NUM_QUERIES 1e6
// #define NUM_QUERIES 1
#define MAX_NUM 10000000
#define NUM_EXPERIMENTS 42

using namespace std;
using namespace std::chrono;
typedef long long ll;
typedef pair<ll, ll> llll;

void clear_cache() {
  ofstream ofs("/proc/sys/vm/drop_caches");
  ofs << "3" << endl;
}

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

vector<int> sorted_array_to_dfs_tree(struct Node *BST, vector<int> &sorted_arr, size_t s) {
  vector<int> result;
  int height = ceil(log2((int)s+1));
  int size = pow(2,height+1)-1;
  result.resize(size);
  for (int i = 0; i < size; i++) result[i] = -1;
  preprocess_dfs_array_rec(BST, 0, height+1, result);
  return result;
}

vector<int> sorted_array_to_bfs_tree(struct Node *BST, vector<int> &sorted_arr, size_t s) {
  vector<int> result;
  int height = ceil(log2((int)s+1));
  int size = pow(2,height+1)-1;
  result.resize(size);
  for (int i = 0; i < size; i++) result[i] = -1;
  preprocess_bfs_array_rec(BST, 0, height+1, result);
  return result;
}

vector<int> sorted_array_to_inorder_tree(struct Node *BST, vector<int> &sorted_arr, size_t s) {
  vector<int> result;
  int height = ceil(log2((int)s+1));
  int size = pow(2,height+1)-1;
  result.resize(size);
  for (int i = 0; i < size; i++) result[i] = -1;
  preprocess_inorder_array_rec(BST, size/2-1, height+1, result);
  return result;
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
  struct Node *BST = sorted_array_to_BST(sorted_arr, 0, s-1);
  // print_array(sorted_arr);
  vector<int> dfs = sorted_array_to_dfs_tree(BST, sorted_arr,s);
  vector<int> bfs = sorted_array_to_bfs_tree(BST, sorted_arr,s);
  vector<int> inorder = sorted_array_to_inorder_tree(BST, sorted_arr,s);
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

llll tester(int(*left)(const int, const int),
                                  int(*right)(const int, const int),
                                  vector<int> &arr,
                                  int height,
                                  int root) {

  PAPI_library_init(PAPI_VER_CURRENT);

  int evts[2] = {PAPI_BR_MSP, PAPI_BR_CN};
  long long values[2];
  
  long long val1 = 0;
  long long val2 = 0;
  for (int ex = 0; ex < NUM_EXPERIMENTS; ex++) {
    clear_cache();

    PAPI_start_counters(evts, 2);

    //run predecessor queries
    for (size_t i = 0; i < NUM_QUERIES; i++) {
      int r = rand() % MAX_NUM;
      int p = tree_predecessor(left, right, arr, r, height, root);
      p++;
    }
    //read & stop
    PAPI_stop_counters(values, 2);
    // cout << values[0] << endl;
    val1 += values[0];
    val2 += values[1];

  }

  return make_pair(val1/NUM_EXPERIMENTS, val2/NUM_EXPERIMENTS);
}

void test_running_time() {
  ofstream lb,bs,dfsf,bfsf,inorderf;
  lb.open("pred_running_time_lb.dat");
  bs.open("pred_running_time_bs.dat");
  dfsf.open("pred_running_time_dfs.dat");
  bfsf.open("pred_running_time_bfs.dat");
  inorderf.open("pred_running_time_inorder.dat");
  //test the running time as a function of input size
  high_resolution_clock::time_point start, end;
  for (double i = 2; i <= 25; i+=1) {
    size_t size = (size_t)pow((double)2, i)-1;
    cout << "Generating test data of size " << size << "..." << endl;
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dis(0, INT_MAX);
    vector<int> data;
    data.resize(size);
    for (size_t j = 0; j < size; j++) {
      data[j] = dis(gen);
    }
    preprocess_sorted_array(data);

    cout << "Running test on lower_bound..." << endl;
    int p,x;
    start = high_resolution_clock::now();
    for (size_t j = 0; j < 1e6; j++) {
      x = dis(gen);
      p = pred_lower_bound(x, data);
    }
    end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start).count();
    lb << size << "\t" << duration << endl;
    cout << "done with test on lower_bound..." << p << endl;

    cout << "Running test on binary search..." << endl;
    start = high_resolution_clock::now();
    for (size_t j = 0; j < 1e6; j++) {
      x = dis(gen);
      p = pred_sorted_array(x, data);
    }
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end - start).count();
    bs << size << "\t" << duration << endl;
    cout << "done with test on binary search..." << p << endl;

    struct Node *BST = sorted_array_to_BST(data, 0, size-1); 

    vector<int> dfs = sorted_array_to_dfs_tree(BST,data,size);
    int dfs_height = ceil(log2((int)dfs.size()));
    cout << "Running test on dfs..." << endl;
    start = high_resolution_clock::now();
    for (size_t j = 0; j < 1e6; j++) {
      x = dis(gen);
      p = tree_predecessor(dfs_left, dfs_right, dfs, x, dfs_height, 0);
    }
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end - start).count();
    dfsf << size << "\t" << duration << endl;
    cout << "done with test on dfs..." << p << endl;
    dfs.clear();

    vector<int> bfs = sorted_array_to_bfs_tree(BST,data,size);
    int bfs_height = ceil(log2((int)bfs.size()));
    cout << "Running test on bfs..." << endl;
    start = high_resolution_clock::now();
    for (size_t j = 0; j < 1e6; j++) {
      x = dis(gen);
      p = tree_predecessor(bfs_left, bfs_right, bfs, x, bfs_height, 0);
    }
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end - start).count();
    bfsf << size << "\t" << duration << endl;
    cout << "done with test on bfs..." << p << endl;
    bfs.clear();

    vector<int> inorder = sorted_array_to_inorder_tree(BST,data,size);
    int inorder_height = ceil(log2((int)inorder.size()));
    int inorder_root = inorder.size()/2-1;
    cout << "Running test on inorder..." << endl;
    start = high_resolution_clock::now();
    for (size_t j = 0; j < 1e6; j++) {
      x = dis(gen);
      p = tree_predecessor(inorder_left, inorder_right, inorder, x, inorder_height, inorder_root);
    }
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end - start).count();
    inorderf << size << "\t" << duration << endl;
    cout << "done with test on bfs..." << p << endl;
    inorder.clear();
  }
}

int main() {
  test_running_time();
  return 0;
  srand (time(NULL));
  //test_bfs();
  //test();
  //return 0;

  vector<pair<int, double > > inorder_res, dfs_res, bfs_res;
  //initialize seed
  for (int ex = 2; ex <= 25; ex++) {
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
    struct Node *BST = sorted_array_to_BST(sorted_arr, 0, s-1);
    vector<int> dfs = sorted_array_to_dfs_tree(BST,sorted_arr,s);
    vector<int> bfs = sorted_array_to_bfs_tree(BST,sorted_arr,s);
    vector<int> inorder = sorted_array_to_inorder_tree(BST,sorted_arr,s);
    // print_array(dfs);
    // print_array(bfs);
    // print_array(inorder);
    int inorder_height = ceil(log2((int)inorder.size()));
    int bfs_height = ceil(log2((int)bfs.size()));
    int dfs_height = ceil(log2((int)dfs.size()));
    int inorder_root = inorder.size()/2-1;

    llll bfs_test = tester(bfs_left, bfs_right, bfs, bfs_height, 0);
    llll dfs_test = tester(dfs_left, dfs_right, dfs, dfs_height, 0);
    llll ino_test = tester(inorder_left, inorder_right, inorder, inorder_height, inorder_root);


    // ratios
    bfs_res.push_back(make_pair(s+1, (double)bfs_test.first/(double)bfs_test.second));
    dfs_res.push_back(make_pair(s+1, (double)dfs_test.first/(double)dfs_test.second));
    inorder_res.push_back(make_pair(s+1, (double)ino_test.first/(double)ino_test.second));
    

    // bfs_res.push_back(make_pair(ex, tester(bfs_left, bfs_right, bfs,
    //                                        bfs_height, 0)));

    
    // dfs_res.push_back(make_pair(ex, tester(dfs_left, dfs_right, dfs,
    //                                        dfs_height,0)));
    
    // inorder_res.push_back(make_pair
    //                       (ex, tester(inorder_left, inorder_right,
    //                                   inorder, inorder_height, inorder_root)));
  }
  cout << "INORDER" << endl;
  print_to_plot(inorder_res, "INO.dat");
  cout << "BFS" << endl;
  print_to_plot(bfs_res, "BFS.dat");
  cout << "DFS" << endl;
  print_to_plot(dfs_res, "DFS.dat");

  return 0;
}
