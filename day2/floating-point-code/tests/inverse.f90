PROGRAM inverse
! a simple program to check how many inverse are not accurate.
! set "prec=8" for testing double precision
  IMPLICIT NONE
  integer, parameter :: prec=4
  REAL(kind=prec) :: X,Y,Z
  INTEGER:: I,J

  i=0
  DO j=1,100,1
     X=REAL(j,prec)
     Y=1.0_prec/X
     Z=Y*X
     IF (Z.ne.1.0_prec) THEN
        WRITE(*,'(A,I4,A,F20.18)') "Inverse for ",j," is not correct=", Z
        i=i+1
     END IF
  END DO
  WRITE(*,*)"found", I , " problematic inverse"
  
END PROGRAM inverse
