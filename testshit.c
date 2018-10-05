#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "2local-qaoa.h"

//*Print out the integer x in binary
int main() {
  int x= 0;
  int n=1;
  int j = 5;
  int i =2;
  printf("Hello, World!\n");
  printf("%d \n",(x>>n));
  printf("%d \n",(x>>n)&1);
  printf("%d \n", (1<<i));

  printf("%d \n", j&(1<<i));


  return 0;
} 


/* int main() {
    printf("Hello, World!");
    return 0;
}
     */