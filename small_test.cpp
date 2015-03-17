#include <iostream>
#include <papi.h>
#include <stdlib.h>
#include <stdio.h>
using namespace std;

int main () {
  int eventset = PAPI_NULL;
  PAPI_library_init(PAPI_VER_CURRENT);
  long long values[1] = {0};
  
  PAPI_create_eventset(&eventset); 

  cout << PAPI_add_event(eventset, PAPI_L1_TCM) << endl;

  cout << PAPI_start(eventset) << endl;
  int i = 0;
  for (int i = 0; i < 100; i++) 
    i *= 3;
  PAPI_read(eventset,values);
  PAPI_stop(eventset,values);
  cout << values[0] << endl;
  return 0;
}
