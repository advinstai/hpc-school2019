#include <stdio.h>
#include <stdlib.h>

#define MAXSIZE 500000000 /* Maximum size of array*/
#define MAXWORKERS 12 /* Maximum amount of worker threads */

int size = MAXSIZE;
int vector[MAXSIZE];
double start_time, end_time; /* start and end times */

void Qsort(int first, int last) {
  int pivot, i_pivot, temp, left, right;
  if (first >= last) return; 

  i_pivot = (first + last) / 2;
  pivot = vector[i_pivot];
  left = first; 
  right = last;
  while (left <= right) {
    if (vector[left] > pivot) { 
       temp = vector[left]; 
       vector[left] = vector[right]; 
       vector[right] = temp;
       if (right == i_pivot) {
        i_pivot = left;
       }
       right--;
    } 
    else { 
      left++;
    }
  }
  // place the pivot in its place (i.e. swap with right element)
  temp = vector[right];
  vector[right] = pivot;
  vector[i_pivot] = temp;
  // sort two sublists in parallel;

    Qsort(first, (right - 1));
    Qsort((right + 1), last);
}

int main(int argc, char *argv[]) {
  int i;

  /* determine size */
  size = (argc > 1) ? atoi(argv[1]) : MAXSIZE;
  if (size <= 0 || size > MAXSIZE)
    size = MAXSIZE;



  /* initialize and print the vector to be sorted */
  for (i = 0; i < size; i++)
  vector[i] = (int) random () % MAXSIZE;


  Qsort(0, (size - 1));

  for (i = 0; i < size - 1; i++)
  if (vector[i] > vector[i + 1]) {
    printf("The resulting vector is not sorted!\n");
  }
  return(0);
}
