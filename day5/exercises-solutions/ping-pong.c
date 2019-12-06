#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv)
{
    int size,rank,mesg1,mesg2;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if (size != 2) {
        if (rank == 0) puts("must have exactly 2 MPI ranks");
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    mesg1 = 42;
    mesg2 = 0;

    if (rank == 0) {
        MPI_Send(&mesg1, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(&mesg2, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Final result is %d\n", mesg2);
    } else { 
        MPI_Recv(&mesg2, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        mesg1 = mesg2 + 99;
        MPI_Send(&mesg1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
