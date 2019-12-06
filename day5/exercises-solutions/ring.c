#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv)
{
    int i,size,rank,buf,mesg,prev,next;
    MPI_Request req;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    prev = rank ? rank - 1 : size - 1;
    next = (rank + 1) % size;

    mesg = rank;
    for (i=0; i < size; ++i) {
        MPI_Irecv(&buf, 1, MPI_INT, prev, 0, MPI_COMM_WORLD, &req);
        MPI_Send(&mesg, 1, MPI_INT, next, 0, MPI_COMM_WORLD);
        MPI_Wait(&req, MPI_STATUS_IGNORE);
        mesg = buf;
        if (rank == 0) printf("result at step %d is: %d\n", i, buf);
    }
    if (rank == 0) printf("final result is: %d\n", buf);
    MPI_Finalize();
    return 0;
}
