#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <algorithm>
#include <vector>
#include <queue>
#include <papi.h>
#include <math.h>

#define SIZE 1024*1024-1
#define NUM_QUERIES 2
#define MAX_NUM 1024 * 1024

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

int bin_search(int(*left)(const int, const int),   //function pointer on how to go left in bst
	       int(*right)(const int, const int),  //function pointer on how to go right in bst
	       int height,                         //current height in the search
	       vector<int> &arr,                   //arr describing the bst
	       int element,                        //the element to search for
	       int root) {                         //current position in tree
  printf("%d %d %d %d %d %d\n", height, element, root, arr[root], left(root,height), right(root,height));
  if (root > SIZE) return -1;
  if (element == arr[root]) return element;
  if (element < arr[root])
    return bin_search(left, right, height-1, arr, element, left(root, height));
  else
    return bin_search(left, right, height-1, arr, element, right(root, height));
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

  while (height) {
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
  return (x+1)*2;
}

int bfs_left(int x, int height) {
  return x*2+1;
}

int dfs_left(int x, int height) {
  return x+1;
}

int dfs_right(int x, int height) {
  return x + pow(2, height - 1);
}

int inorder_left(int x, int height) {
  return x - pow(2, height - 2);
}

int inorder_right(int x, int height) {
  return x + pow(2, height - 2);
}

void print_array(vector<int> a) {
  printf("size %d \n", (int)a.size());
  for (int i = 0; i < a.size(); i++)
    printf("%d ", a[i]);
  puts("");
}

void preprocess_sorted_array(vector<int> &a) {
  sort(a.begin(), a.end());
}

int cur = 0;
void preprocess_dfs_array(struct Node *BST, vector<int> &ret) {
  if (BST == NULL) return;
  ret[cur++] = BST->data;
  preprocess_dfs_array(BST->left, ret);
  preprocess_dfs_array(BST->right, ret);
}

void preprocess_bfs_array(struct Node *BST, vector<int> &bfs_arr) {
  queue<struct Node*> Q;
  int next = 0;
  Q.push(BST);
  while (!Q.empty()) {
    struct Node* t = Q.front(); Q.pop();
    bfs_arr[next++] = t->data;
    if (t->left != NULL)
      Q.push(t->left);
    if (t->right != NULL)
      Q.push(t->right);
  }
}

vector<int> sorted_array_to_dfs_tree(vector<int> &sorted_arr) {
  struct Node *BST = sorted_array_to_BST(sorted_arr, 0, sorted_arr.size());
  vector<int> result;
  result.resize(SIZE);
  preprocess_dfs_array(BST, result);
  return result;
}

vector<int> sorted_array_to_bfs_tree(vector<int> &sorted_arr) {
  struct Node *BST = sorted_array_to_BST(sorted_arr, 0, sorted_arr.size());
  vector<int> result;
  result.resize(SIZE);
  preprocess_bfs_array(BST, result);
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
  //print_array(sorted_arr);
  vector<int> dfs = sorted_array_to_dfs_tree(sorted_arr);
  vector<int> bfs = sorted_array_to_bfs_tree(sorted_arr);
  //vector<int> inorder = sorted_arr;
  int x = rand() % MAX_NUM;
  int height = ceil(log2((int)SIZE));
  int std_pred = pred_lower_bound(x, sorted_arr);
  int dfs_pred = tree_predecessor(dfs_left, dfs_right, dfs, x, height, 0);
  int bfs_pred = tree_predecessor(bfs_left, bfs_right, bfs, x, height, 0);
  //int inorder_pred = tree_predecessor(inorder_left, inorder_right, sorted_arr, x, height, SIZE/2);
  int bin_search = pred_sorted_array(x, sorted_arr);

  printf("%d %d %d %d\n", std_pred, dfs_pred, bfs_pred, bin_search);

}

int main() {
  test();
  return 0;
}

int main3() {
  vector<int> a;
  for (int i = 1; i <= 15; i++) a.push_back(i);
  a.resize(SIZE);
  vector<int> dfs = sorted_array_to_dfs_tree(a);
  print_array(dfs);

  vector<int> bfs = sorted_array_to_bfs_tree(a);
  print_array(bfs);

  printf("tree min: %d\n", dfs[tree_minimum(dfs_left, 8, ceil(log2((int)dfs.size())))-1]);
  printf("tree pre dfs: %d\n", tree_predecessor(dfs_left, dfs_right, dfs, 10, ceil(log2((int)dfs.size())), 0));
  printf("tree pre bfs: %d\n", tree_predecessor(bfs_left, bfs_right, bfs, 10, ceil(log2((int)dfs.size())), 0));
  printf("%d\n", bin_search(dfs_left, dfs_right, ceil(log2((int)dfs.size())), dfs, 1, 0));
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
