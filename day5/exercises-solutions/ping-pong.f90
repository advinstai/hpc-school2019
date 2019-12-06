PROGRAM pingpong
  USE mpi
  IMPLICIT NONE
  INTEGER :: size, rank, ierr, mesg1, mesg2, status(MPI_STATUS_SIZE)

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

  IF (size /= 2) THEN
      IF (rank == 0) PRINT*,'must have exacrtly 2 MPI ranks'
      CALL MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      STOP 'ABORT'
  END IF

  mesg1 = 42
  mesg2 = 0

  IF (rank == 0) THEN
      CALL MPI_SEND(mesg1, 1, MPI_INTEGER, 1, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_RECV(mesg2, 1, MPI_INTEGER, 1, 0, MPI_COMM_WORLD, status, ierr)
      PRINT*,"Final result is", mesg2
  ELSE
      CALL MPI_RECV(mesg2, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, status, ierr)
      mesg1 = mesg2 + 99
      CALL MPI_SEND(mesg1, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
  END IF

  CALL MPI_FINALIZE(ierr)
END PROGRAM pingpong
