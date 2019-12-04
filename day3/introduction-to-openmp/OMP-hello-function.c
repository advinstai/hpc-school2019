#include <stdio.h>
#include <omp.h>

int main() {

    printf(" omp_get_max_threads %d", omp_get_max_threads() );
    printf(" omp_get_thread_limit %d", omp_get_thread_limit() );

    char hn[600];
    int ID = 0;		
    #pragma omp parallel
    {
	ID = omp_get_thread_num();
        gethostname(hn,600);

        printf("hello from hostname %s Thread Number: %d\n",hn, ID);
    }

    printf("Executing with 4 threads");

    #pragma omp parallel num_threads(4)
    {
	ID = omp_get_thread_num();
        gethostname(hn,600);
        printf("\nhello from hostname %s Thread Number: %d\n",hn, ID);
    }

    printf("Executing with 8 threads");
    omp_set_num_threads(8);	
    #pragma omp parallel 
    {
	ID = omp_get_thread_num();
        gethostname(hn,600);
        printf("\nhello from hostname %s Thread Number: %d\n",hn, ID);
    }

    return(0);
}
