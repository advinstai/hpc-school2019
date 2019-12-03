#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "mkl.h"

#define SIZE 12000

double seconds(){
  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
  return sec;
}


int main(int argc, char **argv){
  
  int i, j, k;
  double * A, * B, * C;
  double t_start_mine, t_stop_mine, t_start_dgemm, t_stop_dgemm;
  
  A = (double *) malloc( SIZE * SIZE * sizeof(double) );
  B = (double *) malloc( SIZE * SIZE * sizeof(double) );
  C = (double *) malloc( SIZE * SIZE * sizeof(double) );
  memset( C, 0, SIZE * SIZE * sizeof(double) );

  // Initialization of matrices A and B to random values 
  for( i = 0; i < SIZE * SIZE; i++ ){

    A[i] = (double) rand();
    B[i] = (double) rand();
    
  }
  
  t_start_dgemm = seconds();
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
	      SIZE, SIZE, SIZE, 1.0, A, SIZE, B, SIZE, 0.0, C, SIZE);  
  t_stop_dgemm = seconds();

  //  memset( C, 0, SIZE * SIZE * sizeof(double) );
  t_start_mine = seconds();  
  /*  for( i = 0; i < SIZE; i++ ){
    for( j = 0; j < SIZE; j++ ){
      for( k = 0; k < SIZE; k++ ){
	C[ i * SIZE + j ] += A[ i * SIZE + k ] * B[ k * SIZE + j ];
      }
    }
  }*/
  t_stop_mine = seconds();  

  fprintf( stdout, "\tTime to solution DGEMM of %.3g seconds for SIZE = %d \n", t_stop_dgemm - t_start_dgemm, SIZE );
  fprintf( stdout, "\tTime to solution my bersion of %.3g seconds for SIZE = %d \n", t_stop_mine - t_start_mine, SIZE );
  
  free(A); 
  free(B);
  free(C);

  return 0;
}
