#include <stdio.h>

int main() {

    char hn[600];

    #pragma omp parallel
    {
        gethostname(hn,600);
        printf("hello from hostname %s\n",hn);
    }
    return(0);
}
