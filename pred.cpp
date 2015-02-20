//TODO: create methods for going left and right in dfs & bfs arrays
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <algorithm>
#include <vector>
#include <queue>
#include <papi.h>
#include <math.h>

#define SIZE 15
#define NUM_QUERIES 2
#define MAX_NUM 1024 * 1024 * 1024

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

struct Node* sorted_array_to_BST(vector<int> a, int start, int end) {
  if (start > end) return NULL;
  int mid = (start+end)/2;
  struct Node *root = newNode(a[mid]);
  root->left = sorted_array_to_BST(a, start, mid-1);
  root->right = sorted_array_to_BST(a,mid+1, end);
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

//predecessor query in a bfs layout tree
int pred_bfs_array(int x, vector<int> &bfs_arr) {
  int cur = bfs_arr[0];
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
  printf("%d \n", (int)a.size());
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

int pred_dfs(int x, vector<int> &arr) {
  if (x < arr[0]) return -1;
  int cur_root = 0;
  int height = ceil(log2((int)arr.size()));
  int right = dfs_right(cur_root, height);
  int left = dfs_left(cur_root, height);
  int old_root = 0;
  while (height > 0) {
    printf("%d %d %d %d\n", cur_root, old_root, left, right);
    old_root = cur_root;
    if (arr[right] <= x) cur_root = right;
    else if (arr[left] <= x) cur_root = left;
    else return arr[cur_root];
    --height;
    right = dfs_right(cur_root, height);
    left = dfs_left(cur_root, height);
  }
  return arr[old_root];
}

vector<int> sorted_array_to_dfs_tree(vector<int> sorted_arr) {
  struct Node *BST = sorted_array_to_BST(sorted_arr, 0, sorted_arr.size());
  vector<int> result;
  result.resize(SIZE);
  preprocess_dfs_array(BST, result);
  return result;
}

int main() {

  vector<int> a;
  for (int i = 1; i <= 15; i++) a.push_back(i);
  a.resize(SIZE);
  vector<int> dfs = sorted_array_to_dfs_tree(a);
  print_array(dfs);
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
