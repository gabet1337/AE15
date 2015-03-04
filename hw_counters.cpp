#include "papi.h"
#include <stdlib.h>
#include <stdio.h>

using namespace std;

int main() {

  printf("%d\n",PAPI_num_counters());
  return 0;

}
