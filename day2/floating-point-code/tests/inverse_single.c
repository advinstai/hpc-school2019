

/* a simple program to check how many inverse are not invertible */
#include <stdio.h>

typedef float real;
#define ONE 1.0f

int main(int argc, char **argv)
{

    real x,y,z;
    int i,j;

    i=0;
    for (j=0; j < 100; ++j) {
        x=(real) (j+1);
        y=ONE/x;
        z=y*x;
        if (z != ONE) {
            printf("Inverse for %d is not correct: %12.8f\n",j+1,z);
            ++i;
        }
    }

    printf("Found %d problematic inverse\n", i);
}

