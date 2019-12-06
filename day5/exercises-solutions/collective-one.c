#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv)
{
    int mesg,buf,rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    mesg = rank;
    MPI_Reduce(&mesg,&buf,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    if (rank == 0) printf("final result is: %d\n",buf);
    MPI_Finalize();
    return 0;
}
