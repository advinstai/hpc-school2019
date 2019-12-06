#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
    double tmp;
    int i,j,k;
    int arows,acols;
    double *amatrix;
    int brows,bcols;
    double *bmatrix;
    int crows,ccols;
    double *cmatrix;

    /* read matrix A */
    fscanf(stdin,"%d%d",&arows,&acols);
    amatrix = (double *)malloc(arows*acols*sizeof(double));
    for (i=0; i<arows; ++i) {
        for (j=0; j<acols; ++j) {
            fscanf(stdin,"%lg",&(amatrix[i*arows+j]));
        }
    }
    /* read matrix B */
    fscanf(stdin,"%d%d",&brows,&bcols);
    bmatrix = (double *)malloc(brows*bcols*sizeof(double));
    for (i=0; i<brows; ++i) {
        for (j=0; j<bcols; ++j) {
            fscanf(stdin,"%lg",&(bmatrix[i*brows+j]));
        }
    }

    /* allocate matrix C */
    if ((arows != brows) || (acols != bcols)) {
        printf("Matrices not compatible: %d vs. %d, %d vs. %d\n",
               arows, brows, acols, bcols);
        free(amatrix);
        free(bmatrix);
        return 1;
    }
    crows = arows;
    ccols = acols;
    cmatrix = (double *)malloc(crows*ccols*sizeof(double));

    /* multiplication */
    for (i=0; i<crows; ++i) {
        for (j=0; j<ccols; ++j) {
            tmp = 0.0;
            for (k=0; k < brows; ++k)
                tmp += amatrix[i*arows+k]*bmatrix[j*brows+k];
            cmatrix[i*crows+j] = tmp;
        }
    }

    /* output result */
    for (i=0; i<crows; ++i) {
        for (j=0; j<ccols; ++j) {
            printf("%.15lg ",cmatrix[i*crows+j]);
        }
        printf("\n");
    }
    free(amatrix);
    free(bmatrix);
    free(cmatrix);
    return 0;
}
