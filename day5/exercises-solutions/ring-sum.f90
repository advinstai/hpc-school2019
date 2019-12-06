PROGRAM ring_sum
  USE mpi
  IMPLICIT NONE
  INTEGER :: i, size, rank, buf, mesg, prev, next, ierr, total
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


  mesg = rank
  total = rank
  DO i=2,size
      CALL MPI_IRECV(buf, 1, MPI_INTEGER, prev, 0, MPI_COMM_WORLD, req, ierr)
      CALL MPI_SEND(mesg, 1, MPI_INTEGER, next, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_WAIT(req, status, ierr)
      mesg = buf
      total = total + buf
  END DO

  PRINT*,"final result on rank", rank, "  is:", total

  CALL MPI_FINALIZE(ierr)
END PROGRAM ring_sum
