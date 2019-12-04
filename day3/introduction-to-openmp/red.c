#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#define nt 10
#define N 50000
#define seed 10000

int main(int argc,char **argv){
    double wtime;
    int i,j,max,v[N];

    //generate the numbers of vector with seed
    for (i = 0; i < N; i++) {
        v[i] = i; //rand()%seed;
    }
/*
    //initialize global min
    max=v[0];

    //parallel for, with reduction operation(logN) of result from each thread
    #pragma omp parallel for private(i) shared(v) 
    //reduction(max: max) num_threads(nt)
    for(i =0; i <N; i++)
    {
            if(v[i]>max){
                max=v[i];
            }
    }

    printf("max %d ",max);

*/
	double sum = 0;
	//#pragma omp parallel for reduction(+:sum)
	#pragma omp parallel for  reduction(+:sum)
	for (j = 0; j < N; j++) {
		//#pragma omp atomic
		sum = sum + v[j];
	}

	printf("sum %f ",sum);
    return 0;
}
