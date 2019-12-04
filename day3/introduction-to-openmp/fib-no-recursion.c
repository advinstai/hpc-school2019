#include <stdlib.h>
#include <stdio.h>

#define N 30

int main ( ) {
  unsigned long long fibs[N];

  fibs[0] = fibs[1] = 1;

  int i;
  unsigned long long sum=2;
  #pragma omp parallel for 
  for (i = 2; i < N; i++) {
    fibs[i] = fibs[i - 1] + fibs[i - 2];
    printf("%llu \n", fibs[i]);
  }

  return 0;
}
