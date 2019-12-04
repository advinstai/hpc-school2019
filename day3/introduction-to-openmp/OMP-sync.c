#include <omp.h>
#include <stdio.h>

int main () 
{
  int thid;

  #pragma omp parallel private(thid)
  {
    #pragma omp master
    {
        thid=	omp_get_thread_num();
        printf("master thread only: thread %d \n", thid);        
    }

    thid=	omp_get_thread_num();
    printf("ALL threads: BE CAREFULL! thread  %d \n", thid);        
        
    #pragma omp barrier

    #pragma omp single
    {
        thid=	omp_get_thread_num();
        printf("some thread execute this part (only one): thread  %d \n", thid);        
    }

    #pragma omp barrier

    thid=	omp_get_thread_num();
    printf("after omp barrier! thread %d \n", thid);        
  }

  return 0;
}
