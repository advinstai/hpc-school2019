#include <stdio.h>
int main (int argc, char **argv) {

  int i;
  float x;

  i=0;
  for (x=0.0; x < 1.0; x += 0.01, ++i) {
      printf("x[%d]=%f\n", i, x);
  }
  x=0.0;
  for (i=0; i < 100; x += 0.01, ++i) {
      printf("x[%d]=%f\n", i, x);
  }
  return 0;
}
