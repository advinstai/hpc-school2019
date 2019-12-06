#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

#define TAG 100
#define MPI_PROC_ROOT 0

void read_matrix( double * mat, int n_rows, int n_cols, FILE * fp ){

  int i, j; 

  for( i = 0; i < n_rows; ++i ){
        for( j = 0; j < n_cols; ++j ){
            fscanf( fp, "%lg", &( mat[ i * n_cols + j ] )  );
        }
    }
}

void dist_matrix( double * mat, int n_rows, int n_cols, int rank, int npes, FILE * fp ){

  int count; 
  double * buf_mat = (double *) malloc( n_rows * n_cols * sizeof(double) );  

  if( rank == MPI_PROC_ROOT ){
    read_matrix( mat, n_rows, n_cols, fp );
    for( count = 1; count < npes; count++ ){
      read_matrix( buf_mat, n_rows, n_cols, fp );
      MPI_Send( buf_mat, n_rows * n_cols, MPI_DOUBLE, count, TAG, MPI_COMM_WORLD );
    }
  } 
  else
    MPI_Recv( mat, n_rows * n_cols, MPI_DOUBLE, MPI_PROC_ROOT, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
  
  free( buf_mat );

}


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

    int rank, npes; 
    int loc_size_arows, loc_size_brows;

    FILE * fp = NULL;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &npes) ;

    if( rank == MPI_PROC_ROOT ){ 
      
      fp = fopen( "matrix.dat", "r" );
      
      /* read dimensions matrix A */
      fscanf( fp, "%d%d", &arows, &acols);
      if( arows % npes ){
	fprintf( stderr, "\nDimension no tocmpatible with the number of processes... mod(arows, npes) != 0! \n" );
	MPI_Abort( MPI_COMM_WORLD, 1 );
      }
    }      

    MPI_Bcast( &arows, 1, MPI_INT, MPI_PROC_ROOT, MPI_COMM_WORLD );
    MPI_Bcast( &acols, 1, MPI_INT, MPI_PROC_ROOT, MPI_COMM_WORLD );

    /* read matrix and distribute matrix A */
    loc_size_arows = arows / npes;
    amatrix = (double *) malloc( loc_size_arows * acols * sizeof(double) );
    dist_matrix( amatrix, loc_size_arows, acols, rank, npes, fp );

    if( rank == MPI_PROC_ROOT ){ 

      /* read dimensions matrix B */
      fscanf(fp,"%d%d", &brows, &bcols);
      if( brows % npes ){
	fprintf( stderr, "\nDimension no tocmpatible with the number of processes... mod( brows, npes ) != 0! \n" );
	MPI_Abort( MPI_COMM_WORLD, 1 );
      }

    }      

    MPI_Bcast( &brows, 1, MPI_INT, MPI_PROC_ROOT, MPI_COMM_WORLD );
    MPI_Bcast( &bcols, 1, MPI_INT, MPI_PROC_ROOT, MPI_COMM_WORLD );

    /* read matrix and distribute matrix B */
    loc_size_brows = brows / npes;
    bmatrix = (double *) malloc( loc_size_brows * bcols * sizeof(double) );
    dist_matrix( bmatrix, loc_size_brows, bcols, rank, npes, fp );

    /* allocate matrix C */
    if ((arows != brows) || (acols != bcols)) {
        printf("Matrices not compatible: %d vs. %d, %d vs. %d\n",
               arows, brows, acols, bcols);
        free(amatrix);
        free(bmatrix);
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    crows = loc_size_arows;
    ccols = acols;
    cmatrix = (double *) malloc( crows * ccols * sizeof(double) );
    memset( cmatrix, 0, crows * ccols * sizeof(double) );

    /* multiplication */
    /* Implement here the parallel matrix matmul */

    free(amatrix);
    free(bmatrix);
    free(cmatrix);
    if( fp ) fclose( fp );
    
    MPI_Finalize();
    return 0;
}
