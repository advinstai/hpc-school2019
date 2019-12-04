#include <stdio.h>
#include <stdlib.h>

double *b;

void print(int n) {
  int i;
  for(i=0; i<n; i++)
    printf("%.2lf ", b[i]);
  printf("\n");
}

void metDARE(int n, double alpha, double beta, double T, int k) {
  int i, j;
  
  b = (double *) malloc(sizeof(double)*n);
  b[n-1] = T;

  #pragma omp parallel for
  for(i=0; i<n-1; i++)
    b[i] = 0.0;

  for(j=0; j<k; j++) {
    for(i=0; i<n-1; i++) {
      b[i] = 0.5*(alpha*b[i] + beta*b[i+1]);
    }
    b[n-1] = b[n-1] + 0.0001; 
  }
  print(n);
  free(b);
}

int main() {
  int n, k;
  double alpha, beta, T;
  scanf("%d %lf %lf %lf %d", &n, &alpha, &beta, &T, &k);
  metDARE(n, alpha, beta, T, k);
  
  return 0;
}
