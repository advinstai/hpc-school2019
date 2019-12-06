PROGRAM hello
  USE mpi
  IMPLICIT NONE
  INTEGER :: i, size, rank, ierr, status(MPI_STATUS_SIZE)

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

  IF (rank == 0) THEN
      DO i=1,size
          IF (i > 1) THEN
              CALL MPI_RECV(rank,1,MPI_INTEGER,i-1,0,MPI_COMM_WORLD,status,ierr)
          END IF
          PRINT*,"Hello, World! I am rank ",rank," out of a total of",size," tasks"
      END DO
  ELSE
      CALL MPI_SEND(rank,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)
  END IF
  CALL MPI_FINALIZE(ierr)
END PROGRAM hello
