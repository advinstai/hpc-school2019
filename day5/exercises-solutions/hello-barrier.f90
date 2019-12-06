PROGRAM hello_barrier
  USE mpi
  IMPLICIT NONE
  INTEGER :: i, size, rank, ierr

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

  DO i=1,size
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      IF (rank == i-1) THEN
          PRINT*,"Hello, World! I am rank ",rank," out of a total of",size," tasks"
      END IF
  END DO
  CALL MPI_FINALIZE(ierr)
END PROGRAM hello_barrier
