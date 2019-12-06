PROGRAM hello
  USE mpi
  IMPLICIT NONE
  INTEGER :: size, rank, ierr

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
  
  PRINT*,"Hello, World! I am rank ",rank," out of a total of",size," tasks"

  CALL MPI_FINALIZE(ierr)
END PROGRAM hello
