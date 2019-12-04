#include <stdio.h>
#include <stdlib.h>

int main(){

  int nn, i;
  nn = 2000000000;
  double* a;
  double* b;
  double* c;

  a = (double*) malloc(nn*sizeof(double));
  b = (double*) malloc(nn*sizeof(double));
  c = (double*) malloc(nn*sizeof(double));

  #pragma omp parallel for
  for( i=0; i<nn; i++) {
    a[i] = i+i;
    b[i] = i-i;
    c[i] = i*i;
  }

  #pragma omp parallel for
  for( i=0; i<nn; i++) {
    a[i] = b[i] + c[i];
  }

  #pragma omp parallel for
  for( i=0; i<nn; i++) {
    b[i] = c[i] + a[i];
  }

  #pragma omp parallel for
  for( i=0; i<nn; i++) {
    c[i] = b[i] + a[i];
  }

}
