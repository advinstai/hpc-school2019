PROGRAM ring
  USE mpi
  IMPLICIT NONE
  INTEGER :: i, size, rank, buf, mesg, prev, next, ierr, dims, periods
  INTEGER :: req, status(MPI_STATUS_SIZE), comm

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
  dims = size
  periods = 1

  CALL MPI_CART_CREATE(MPI_COMM_WORLD, 1, dims, periods, 1, comm, ierr)
  CALL MPI_COMM_SIZE(comm,size,ierr)
  CALL MPI_COMM_RANK(comm,rank,ierr)

  prev = 0
  next = 0
  CALL MPI_CART_SHIFT(comm,0,1,prev,next,ierr)

  mesg = rank
  DO i=1,size
      CALL MPI_IRECV(buf, 1, MPI_INTEGER, prev, 0, comm, req, ierr)
      CALL MPI_SEND(mesg, 1, MPI_INTEGER, next, 0, comm, ierr)
      CALL MPI_WAIT(req, status, ierr)
      mesg = buf
      IF (rank == 0) PRINT*,"result at step", i, "  is:", buf
  END DO

  IF (rank == 0) PRINT*,"final result is:", buf

  CALL MPI_FINALIZE(ierr)
END PROGRAM ring
