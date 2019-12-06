PROGRAM collective_one
  USE mpi
  IMPLICIT NONE
  INTEGER :: rank, buf, mesg, ierr

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

  mesg = rank
  CALL MPI_REDUCE(mesg, buf, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  IF (rank == 0) PRINT*,"final result is:", buf

  CALL MPI_FINALIZE(ierr)
END PROGRAM collective_one
