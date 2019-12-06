PROGRAM ring_once
  USE mpi
  IMPLICIT NONE
  INTEGER :: size, rank, buf, mesg, prev, next, ierr
  INTEGER :: req, status(MPI_STATUS_SIZE)

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

  IF (rank == 0) THEN
      prev = size-1
  ELSE
      prev = rank-1
  END IF
  next = MOD(rank+1,size)


  mesg = 0
  CALL MPI_IRECV(buf, 1, MPI_INTEGER, prev, 0, MPI_COMM_WORLD, req, ierr)

  IF (rank == 0) THEN
      CALL MPI_SEND(mesg, 1, MPI_INTEGER, next, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_WAIT(req, status, ierr)
  ELSE
      CALL MPI_WAIT(req, status, ierr)
      mesg = buf + rank
      CALL MPI_SEND(mesg, 1, MPI_INTEGER, next, 0, MPI_COMM_WORLD, ierr)
  END IF
  IF (rank == 0) PRINT*,"final result is:", buf

  CALL MPI_FINALIZE(ierr)
END PROGRAM ring_once
