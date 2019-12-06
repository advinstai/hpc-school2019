#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv)
{
    int i,size,rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if (rank == 0) {
        for (i=0; i < size; ++i) {
	    if (i > 0)
                MPI_Recv(&rank, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Hello, World! I am rank %d "
                   "out of a total of %d tasks\n", rank, size);
        }
    } else {
      MPI_Send(&rank,1,MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}
