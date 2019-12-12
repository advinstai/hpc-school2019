#include "stdio.h"
#include "stdlib.h"
#include "timer.h"

////#define SIZE 1024
#define SIZE 3072
//#define SIZE 6144

void multiply(float *restrict A, float *restrict B, float *restrict C) 
{
  int i;
  int row, col, row_len=SIZE, col_len=SIZE;


  #pragma acc kernels copyin(A[:SIZE*SIZE],B[:SIZE*SIZE]), copyout(C[:SIZE*SIZE])
  { 
    for (row=0; row<row_len; row++) {
      for (col=0; col<col_len; col++) {
	for (i=0; i<SIZE; i++)
	  C[row*row_len+col] += A[row*row_len+i] * B[col+i*row_len];
      }
    }
  }
}

int main(int argc, char **argv) 
{
  int i;
  float *restrict A;
  float *restrict B;
  float *restrict C;
  
  A = malloc(SIZE * SIZE * sizeof(float));
  B = malloc(SIZE * SIZE * sizeof(float));
  C = malloc(SIZE * SIZE * sizeof(float));

  for (i=0; i<SIZE*SIZE; i++) {
    A[i] = (float)i;
    B[i] = (float)SIZE * (float)SIZE - (float)i - 1.0;
  }

  StartTimer();
  
  multiply(A, B, C); // perform matrix multiply


  double runtime = GetTimer();
  printf("\nElapsed time: %f s\n\n", runtime / 1000.f);

  return 0;
}
