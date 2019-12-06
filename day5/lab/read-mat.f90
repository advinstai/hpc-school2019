PROGRAM matmul

  use mpi

  IMPLICIT NONE

  REAL(kind=8),ALLOCATABLE :: amatrix(:,:), bmatrix(:,:), cmatrix(:,:), buf_mat(:,:)
  REAL(kind=8) :: tmp
  INTEGER :: arows,acols,brows,bcols,crows,ccols
  INTEGER :: i,j,k
  INTEGER :: loc_size_arows, loc_size_brows, loc_size_crows

  INTEGER :: ierror
  INTEGER :: rank, npes, count

  CALL MPI_INIT(ierror)
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierror )
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, npes, ierror )

  ! read dimensions matrix A
  IF( rank == 0 ) THEN 

     open( unit = 20, file = "matrix.dat", status = 'old', form = 'formatted' )  

      READ (20,*) arows,acols

      IF( MOD(arows, npes) /= 0 ) THEN 
         PRINT *, "Dimension no tocmpatible with the number of processes... mod(arows, npes) != 0!" 
         CALL MPI_ABORT( MPI_COMM_WORLD, 1, ierror )
      ENDIF
   ENDIF

   CALL MPI_BCAST( arows, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror )
   CALL MPI_BCAST( acols, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror )

   loc_size_arows = arows / npes
   ALLOCATE( amatrix( loc_size_arows, acols) )

   ! distribute matrix A
   IF ( rank == 0 ) THEN 
      
      ALLOCATE( buf_mat( loc_size_arows, acols) )
      DO i=1,loc_size_arows
         READ (20, *), (amatrix(j,i),j=1,acols)     
      END DO

      do count = 1, npes - 1
         DO i=1,loc_size_arows
            READ (20, *), (buf_mat(j,i),j=1,acols)     
         END DO
         CALL MPI_SEND( buf_mat, loc_size_arows * acols, MPI_DOUBLE_PRECISION, count, 100, MPI_COMM_WORLD, ierror )
      end do
      DEALLOCATE( buf_mat )
      
   ELSE
      CALL MPI_RECV( amatrix, loc_size_arows * acols, MPI_DOUBLE_PRECISION, 0, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror )
   ENDIF

  ! read dimensions matrix B
  IF( rank == 0 ) THEN 

      READ (20,*) brows,bcols

      IF( MOD(brows, npes) /= 0 ) THEN 
         PRINT *, "Dimension no tocmpatible with the number of processes... mod(arows, npes) != 0!" 
         CALL MPI_ABORT( MPI_COMM_WORLD, 1, ierror )
      ENDIF
   ENDIF

   CALL MPI_BCAST( brows, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror )
   CALL MPI_BCAST( bcols, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror )

   loc_size_brows = brows / npes
   ALLOCATE( bmatrix( loc_size_brows, bcols) )

   ! distribute matrix B
   IF ( rank == 0 ) THEN 
      
      ALLOCATE( buf_mat( loc_size_brows, bcols) )
      DO i=1,loc_size_brows
         READ (20, *), (bmatrix(j,i),j=1,bcols)
      END DO

      do count = 1, npes - 1
         DO i=1,loc_size_brows
            READ (20, *), (buf_mat(j,i),j=1,bcols)     
         END DO
         CALL MPI_SEND( buf_mat, loc_size_brows * bcols, MPI_DOUBLE_PRECISION, count, 100, MPI_COMM_WORLD, ierror )
      end do
      DEALLOCATE( buf_mat )
      
   ELSE
      CALL MPI_RECV( bmatrix, loc_size_brows * bcols, MPI_DOUBLE_PRECISION, 0, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror )
   ENDIF

  crows = loc_size_arows
  ccols = acols
  ALLOCATE(cmatrix(crows,ccols))

  ! IMPLEMENT HERE!!

  CALL MPI_FINALIZE( ierror )

END PROGRAM matmul



