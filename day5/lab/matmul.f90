PROGRAM matmul
  REAL(kind=8),ALLOCATABLE :: amatrix(:,:), bmatrix(:,:), cmatrix(:,:)
  REAL(kind=8) :: tmp
  INTEGER :: arows,acols,brows,bcols,crows,ccols
  INTEGER :: i,j,k

  ! read matrix A
  READ*,arows,acols
  ALLOCATE(amatrix(arows,acols))
  DO i=1,arows
      READ*,(amatrix(j,i),j=1,acols)
  END DO
  ! read matrix B
  READ*,brows,bcols
  ALLOCATE(bmatrix(brows,bcols))
  DO i=1,brows
      READ*,(bmatrix(j,i),j=1,bcols)
  END DO
  ! allocate matrix C
  IF ((arows /= brows) .AND. (acols /= bcols)) THEN
      PRINT*, 'Matrices are not compatible. ',arows,' vs ',brows,', ', acols, ' vs ',bcols

      DEALLOCATE(amatrix,bmatrix)
      STOP 'input data error'
  END IF
  crows = arows
  ccols = acols
  ALLOCATE(cmatrix(arows,acols))

  DO i=1,crows
      DO j=1,ccols
          tmp = 0.0_8
          DO k=1,brows
              tmp = tmp + amatrix(k,i)*bmatrix(k,j)
          END DO
          cmatrix(j,i) = tmp
      END DO
  END DO

  ! output result
  DO i=1,crows
      WRITE(*,*) (cmatrix(j,i),j=1,ccols)
  END DO


END PROGRAM matmul

