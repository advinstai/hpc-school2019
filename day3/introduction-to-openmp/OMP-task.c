#include <stdio.h>
int main() {
  #pragma omp parallel​
  {
    #pragma omp single​
    {
      #pragma omp task
      printf("hello \n");

      #pragma omp task
      printf("again \n");
    }
  }
  return(0);
}
