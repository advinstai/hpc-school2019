#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
    int i,size,rank,buf,mesg,prev,next;
    MPI_Request req;
    MPI_Comm comm;
    int dims,periods=1;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    dims = size;

    MPI_Cart_create(MPI_COMM_WORLD, 1, &dims, &periods, 1, &comm);
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    prev = 0;
    next = 0;
    MPI_Cart_shift(comm,0,1,&prev,&next);

    mesg = rank;
    for (i=0; i < size; ++i) {
        MPI_Irecv(&buf, 1, MPI_INT, prev, 0, comm, &req);
        MPI_Send(&mesg, 1, MPI_INT, next, 0, comm);
        MPI_Wait(&req, MPI_STATUS_IGNORE);
        mesg = buf;
        if (rank == 0) printf("result at step %d is: %d\n", i, buf);
    }
    if (rank == 0) printf("final result is: %d\n", buf);
    MPI_Finalize();
    return 0;
}
