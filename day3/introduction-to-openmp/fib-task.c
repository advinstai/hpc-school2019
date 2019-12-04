#include <stdlib.h>
#include <stdio.h>

unsigned long long fib(int n) {

  unsigned long long i, j, res;

  if (n < 2) return n;
  
  #pragma omp task shared(i)
  i = fib(n-1);
  #pragma omp task shared(j)
  j = fib(n-2);

  #pragma omp taskwait
  return (i+j);
}

int main(){
  int n = 80;
  int c = 1;
  #pragma omp parallel
  {
    #pragma omp single
    {
      for (c = 1; c <= n; c++) {
        printf("%llu\n", fib(c));
      }
    }
  }
}
