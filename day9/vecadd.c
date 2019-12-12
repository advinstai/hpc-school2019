#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
int main( int argc, char* argv[] )
{
 
    // Size of vectors
    int n = 100000;
 
    // Input vectors
    double *restrict a;
    double *restrict b;
    // Output vector
    double *restrict c;
 
    // Size, in bytes, of each vector
    size_t bytes = n*sizeof(double);
 
    // Allocate memory for each vector
    a = (double*)malloc(bytes);
    b = (double*)malloc(bytes);
    c = (double*)malloc(bytes);
 
    // Initialize content of input vectors, vector a[i] = sin(i)^2 vector b[i] = cos(i)^2
    int i;
    for(i=0; i<n; i++) {
        a[i] = sin(i)*sin(i);
        //b[i] = cos(i)*cos(i);
	b[i] = tan(i)*tan(i);
        c[i] = tan(i)*tan(i);
    }  

    int nn=20000;
    int j=0;
    for(j=0; j<nn; j++) { 
    // sum component wise and save result into vector c
	#pragma acc kernels copyin(a[0:n],b[0:n]), copyout(c[0:n])
	//#pragma acc parallel  
    	{
	    for(i=0; i<n; i++) {
        	c[i] = a[i] + b[i];
	    }
    	}
    }

    // Sum up vector c and print result divided by n, this should equal 1 within error
    //double sum = 0.0;
    //for(i=0; i<n; i++) {
    //    sum += c[i];
    //}
    //sum = sum/n;
    //printf("final result: %f\n", sum);
 
    // Release memory
    free(a);
    free(b);
    free(c);
 
    return 0;
}
