#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>

int main() {

  double *a;
  double *b;
  double randV;
  double lastx;
  int i =0; 
  int msize = 0;

  //msize = 987456213; //16000;
  msize = 1000; //16000;
  randV = 0;

  srand(time(NULL));
  randV=rand();
  randV=randV;

  a = (double*) malloc(msize * sizeof(double));
  for ( i = 0; i < msize; i++) {
    a[i] = i/randV;
  }

  randV=rand();
  randV=randV;
  printf("randV %f\n", randV);

  b = (double*) malloc(msize * sizeof(double));

  int totthreads;
  
  totthreads=omp_get_num_threads();
  totthreads=11;

  int id=0;
  int stepPThread;    
  int stepLastThread;    
  int begin;
  int end;

  #pragma omp parallel shared(totthreads) private(id, stepPThread, stepLastThread, begin, end)
  {
    id=omp_get_thread_num();
    stepPThread=msize/totthreads;    
    stepLastThread=msize % totthreads;    
    begin = id*stepPThread;
    end = (begin + stepPThread)-1;

    int cmp = totthreads-1;

    if (id == cmp) {
      end = end + stepLastThread;
      //printf("\n id %d cmp %d \n", id, cmp);       
    }
    //else
      //printf("\n id %d cmp %d \n", id, cmp);

    printf("totthreads %d id %d stepPThread %d stepLastThread %d begin %d end %d \n", totthreads, id, stepPThread, stepLastThread, begin, end);

    for ( i = begin; i < end; i++) {
      b[i] = i/randV;
    }
  }

  #pragma omp parallel 
  {
    #pragma omp for
    for ( i = 0; i < msize; i++) {
      b[i] = i/randV;
    }
  }

  #pragma omp parallel for
  for ( i = 0; i < msize; i++) {
    a[i] = a[i]+b[i];
  }

  lastx = 0;
  for ( i = 0; i < msize; i++) {
      lastx = a[i]+b[i] ;
  }

  printf("lastx SERIAL %f\n", lastx);

  lastx = 0;
  #pragma omp parallel for lastprivate(lastx)
  for ( i = 0; i < msize; i++) {
      lastx = a[i]+b[i] ;
  }

  printf("lastx lastprivate %f\n", lastx);

  lastx = 0;
  #pragma omp parallel for
  for ( i = 0; i < msize; i++) {
      lastx =  a[i]+b[i] ;
  }

  printf("lastx omp parallel for %f\n", lastx);
 
  return(0);
}
