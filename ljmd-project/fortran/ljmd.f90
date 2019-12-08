! 
! simple lennard-jones potential MD code with velocity verlet.
! units: Length=Angstrom, Mass=amu, Energy=kcal
!
! baseline f95 version.
!

MODULE kinds
  IMPLICIT NONE
  INTEGER, PARAMETER :: dbl = selected_real_kind(14,200)  ! double precision floating point
  INTEGER, PARAMETER :: sgl = selected_real_kind(6,30)    ! single precision floating point
  INTEGER, PARAMETER :: sln = 200                         ! length of I/O input line
  PRIVATE
  PUBLIC :: sgl, dbl, sln
END MODULE kinds

MODULE physconst
  USE kinds
  IMPLICIT NONE
  REAL(kind=dbl), PARAMETER :: kboltz =    0.0019872067_dbl   ! boltzman constant in kcal/mol/K
  REAL(kind=dbl), PARAMETER :: mvsq2e = 2390.05736153349_dbl  ! m*v^2 in kcal/mol
  PRIVATE
  PUBLIC :: kboltz, mvsq2e
END MODULE physconst

! module to hold the complete system information 
MODULE mdsys
  USE kinds
  IMPLICIT NONE
  INTEGER :: natoms,nfi,nsteps
  REAL(kind=dbl) dt, mass, epsilon, sigma, box, rcut
  REAL(kind=dbl) ekin, epot, temp
  REAL(kind=dbl), POINTER, DIMENSION (:) :: rx, ry, rz
  REAL(kind=dbl), POINTER, DIMENSION (:) :: vx, vy, vz
  REAL(kind=dbl), POINTER, DIMENSION (:) :: fx, fy, fz
END MODULE mdsys

MODULE utils
  USE kinds
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: azzero, pbc

CONTAINS
! helper function: zero out an array 
  SUBROUTINE azzero(d, n)
    REAL(kind=dbl), DIMENSION(:), INTENT(INOUT) :: d
    INTEGER, INTENT(IN) :: n
    INTEGER :: i

    DO i=1, n
       d(i) = 0.0_dbl
    END DO
  END SUBROUTINE azzero
    
! helper function: apply minimum image convention 
  FUNCTION pbc(x, box)
    REAL(kind=dbl), INTENT(IN)  :: x, box
    REAL(kind=dbl) :: pbc

    pbc = x - box*(ANINT(x/box))
  END FUNCTION pbc
END MODULE utils

MODULE io
  USE kinds
  IMPLICIT NONE
  PRIVATE 
  INTEGER, PARAMETER :: stdin=5, stdout=6, log=30, xyz=31
  PUBLIC :: ioopen, ioclose, output, stdin, stdout, getline

CONTAINS
  SUBROUTINE getline(chan, line)
    INTEGER, INTENT(IN) :: chan
    CHARACTER(LEN=sln), INTENT(OUT) :: line
    INTEGER :: idx, i

    READ(CHAN, '(A)') line
    ! delete comment
    idx=INDEX(line,'#')
    IF (idx > 0) THEN
       DO i=idx,sln
          line(i:i) = ' '
       END DO
    END IF
  END SUBROUTINE getline

  SUBROUTINE ioopen(logname, xyzname)
    CHARACTER(LEN=sln) :: logname, xyzname
    OPEN(UNIT=log, FILE=TRIM(logname), STATUS='UNKNOWN', FORM='FORMATTED')
    OPEN(UNIT=xyz, FILE=TRIM(xyzname), STATUS='UNKNOWN', FORM='FORMATTED')
  END SUBROUTINE ioopen
  
  SUBROUTINE ioclose
    CLOSE(UNIT=log)
    CLOSE(UNIT=xyz)
  END SUBROUTINE ioclose
  
  ! append data to output.
  SUBROUTINE output
    USE mdsys
    IMPLICIT NONE
    INTEGER :: i
    WRITE(log, '(I8,1X,F20.8,1X,F20.8,1X,F20.8,1X,F20.8)') &
         nfi, temp, ekin, epot, ekin+epot
    WRITE(stdout, '(I8,1X,F20.8,1X,F20.8,1X,F20.8,1X,F20.8)') &
         nfi, temp, ekin, epot, ekin+epot
    WRITE(xyz, '(I8)') natoms
    WRITE(xyz, '(A,I8,1X,A,F20.8)') 'nfi=', nfi, 'etot=', ekin+epot
    DO i=1, natoms
       WRITE(xyz, '(A,1X,F20.8,1X,F20.8,1X,F20.8)') &
            'Ar ', rx(i), ry(i), rz(i)
    END DO
  END SUBROUTINE output
END MODULE io

! compute kinetic energy
SUBROUTINE getekin
  USE kinds
  USE mdsys, ONLY: natoms, mass, temp, ekin, vx, vy, vz
  USE physconst
  IMPLICIT NONE

  INTEGER :: i

  ekin = 0.0_dbl
  DO i=1, natoms
     ekin = ekin + 0.5_dbl * mvsq2e * mass * (vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i))
  END DO

  temp = 2.0_dbl * ekin/(3.0_dbl*DBLE(natoms-1))/kboltz
END SUBROUTINE getekin

! compute forces 
SUBROUTINE force
  USE kinds
  USE utils
  USE mdsys
  IMPLICIT NONE

  REAL(kind=dbl) :: r, ffac, dx, dy, dz
  INTEGER :: i, j
  
  epot=0.0_dbl
  CALL azzero(fx,natoms)
  CALL azzero(fy,natoms)
  CALL azzero(fz,natoms)

  DO i=1, natoms
     DO j=1, natoms
        ! particles have no interactions with themselves 
        IF (i==j) CYCLE
            
        ! get distance between particle i and j 
        !        delta = delta - box*(ANINT(delta/box))
        dx=pbc(rx(i) - rx(j), box)
        dy=pbc(ry(i) - ry(j), box)
        dz=pbc(rz(i) - rz(j), box)
        r = SQRT(dx*dx + dy*dy + dz*dz)
      
        ! compute force and energy if within cutoff */
        IF (r < rcut) THEN
           ffac = -4.0_dbl*epsilon*(-12.0_dbl*((sigma/r)**12)/r   &
                +6.0_dbl*(sigma/r)**6/r)
                
           epot = epot + 0.5_dbl*4.0_dbl*epsilon*((sigma/r)**12 &
                -(sigma/r)**6.0)

           fx(i) = fx(i) + dx/r*ffac
           fy(i) = fy(i) + dy/r*ffac
           fz(i) = fz(i) + dz/r*ffac
        END IF
     END DO
  END DO
END SUBROUTINE force


! velocity verlet
SUBROUTINE velverlet
  USE kinds
  USE mdsys
  USE physconst
  IMPLICIT NONE

  INTEGER :: i


  ! first part: propagate velocities by half and positions by full step
  DO i=1, natoms
     vx(i) = vx(i) + 0.5_dbl * dt / mvsq2e * fx(i) / mass
     vy(i) = vy(i) + 0.5_dbl * dt / mvsq2e * fy(i) / mass
     vz(i) = vz(i) + 0.5_dbl * dt / mvsq2e * fz(i) / mass
     rx(i) = rx(i) + dt*vx(i)
     ry(i) = ry(i) + dt*vy(i)
     rz(i) = rz(i) + dt*vz(i)
  END DO

  ! compute forces and potential energy 
  CALL force

  ! second part: propagate velocities by another half step */
  DO i=1, natoms
     vx(i) = vx(i) + 0.5_dbl * dt / mvsq2e * fx(i) / mass
     vy(i) = vy(i) + 0.5_dbl * dt / mvsq2e * fy(i) / mass
     vz(i) = vz(i) + 0.5_dbl * dt / mvsq2e * fz(i) / mass
  END DO
END SUBROUTINE velverlet


PROGRAM LJMD
  USE kinds
  USE io
  USE utils
  USE mdsys
  IMPLICIT NONE
  
  INTEGER :: nprint, i
  CHARACTER(len=sln) :: restfile, trajfile, ergfile

  READ(stdin,*) natoms
  READ(stdin,*) mass
  READ(stdin,*) epsilon
  READ(stdin,*) sigma
  READ(stdin,*) rcut
  READ(stdin,*) box
  CALL getline(stdin,restfile)
  CALL getline(stdin,trajfile)
  CALL getline(stdin,ergfile)
  READ(stdin,*) nsteps
  READ(stdin,*) dt
  READ(stdin,*) nprint

  ! allocate storage for simulation data.
  ALLOCATE(rx(natoms),ry(natoms),rz(natoms),&
       vx(natoms),vy(natoms),vz(natoms), &
       fx(natoms),fy(natoms),fz(natoms))

  ! read restart 
  OPEN(UNIT=33, FILE=restfile, FORM='FORMATTED', STATUS='OLD')
  DO i=1,natoms
     READ(33,*) rx(i), ry(i), rz(i)
  END DO
  DO i=1,natoms
     READ(33,*) vx(i), vy(i), vz(i)
  END DO
  CLOSE(33)

  CALL azzero(fx,natoms)
  CALL azzero(fy,natoms)
  CALL azzero(fz,natoms)

  ! initialize forces and energies
  nfi=0
  CALL force
  CALL getekin
    
  CALL ioopen(ergfile, trajfile)

  WRITE(stdout, *) 'Starting simulation with ', natoms, ' atoms for', nsteps, ' steps'
  WRITE(stdout, *) '    NFI           TEMP                 EKIN                  EPOT&
       &                ETOT'
  CALL output

  ! main MD loop 
  DO nfi=1, nsteps
     ! write output, if requested
     IF (mod(nfi,nprint) == 0) THEN
        CALL output
     END IF

        ! propagate system and recompute energies
        CALL velverlet
        CALL getekin
     END DO

     ! clean up: close files, free memory
     WRITE(stdout,'(A)') 'Simulation Done.'
     CALL ioclose

     DEALLOCATE(rx,ry,rz,vx,vy,vz,fx,fy,fz)
END PROGRAM LJMD
