#include <stdio.h>
#include <stdlib.h>

int main(){
    int* input;    
    int snum, i, sum;

    sum = 0;
    snum=40000000;
    
    input = (int*) malloc (sizeof(int)*snum);

    for(i=0;i<snum;i++) {
        input[i] = i+1;
    }

    #pragma omp parallel for
    for(i=0;i<snum;i++)
    {
        int* tmpsum = input+i;
                
        //#pragma omp critical 
        //{
        //    sum += *tmpsum;
        //}
        #pragma omp atomic
        sum += *tmpsum;

    }
    printf("sum %d", sum);
}
