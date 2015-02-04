//TODO: create methods for going left and right in dfs & bfs arrays
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <algorithm>
#include <vector>
#include <queue>

#define SIZE 16
#define NUM_QUERIES 2
#define MAX_NUM 1024 * 1024

using namespace std;

vector<int> sorted_arr;
vector<int> dfs_arr;
vector<int> bfs_arr;

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
int pred_sorted_array(int x) {

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
int pred_lower_bound(int x) {
  vector<int>::iterator it = lower_bound(sorted_arr.begin(), sorted_arr.end(), x);
  if (*it == x) return x;
  if (it == sorted_arr.begin()) return -1;
  return *(it-1);
}

//predecessor query in a bfs layout tree
int pred_bfs_array(int x) {
  
}

void print_array(vector<int> a) {
  printf("%d \n", a.size());
  for (int i = 0; i < a.size(); i++)
    printf("%d ", a[i]);
  puts("");
}

void preprocess_sorted_array() {
  sort(sorted_arr.begin(),sorted_arr.end());
  print_array(sorted_arr);
}

int cur = 0;
void preprocess_dfs_array(struct Node *BST) {
  if (BST == NULL) return;
  dfs_arr[cur++] = BST->data;
  preprocess_dfs_array(BST->left);
  preprocess_dfs_array(BST->right);
}

void preprocess_bfs_array(struct Node *BST) {
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

int main () {
  //initialize seed
  srand (time(NULL));
  //  fill the sorted_array with random numbers
  sorted_arr.resize(SIZE);
  bfs_arr.resize(SIZE);
  dfs_arr.resize(SIZE);

  for (int i = 0; i < SIZE; i++) {
    int r = rand() % MAX_NUM;
    sorted_arr[i] = r;
  }
  
  preprocess_sorted_array();
  struct Node *BST = sorted_array_to_BST(sorted_arr, 0, sorted_arr.size());

  preprocess_dfs_array(BST);
  print_array(dfs_arr); 
  preprocess_bfs_array(BST);
  print_array(bfs_arr);

  //make random queries using binary search in sorted sorted_array
  for (int i = 0; i < NUM_QUERIES; i++) {
    int x = rand() % MAX_NUM; //rand below SIZE
    printf("pred of %d is %d\n", x, pred_sorted_array(x));
  }

  //make random queries using binary search in sorted sorted_array
  for (int i = 0; i < NUM_QUERIES; i++) {
    int x = rand() % MAX_NUM; //rand below SIZE
    printf("pred of %d is %d\n", x, pred_lower_bound(x));
  }

  return 0;
}
