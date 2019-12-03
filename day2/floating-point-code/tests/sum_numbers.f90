   program sum_numbers
! 
! a program to show different results in summing FP numbers
! REMEMBER: FP arithmetic is NOT commutative 
!
! S.C. 28/01/07 for the ICTP activity smr1830
!

   real :: directly=0.0  
   real :: reversly=0.0
   integer :: i_dir=0
   integer :: i_rev=0
   integer :: nsteps ,i,ilarge
   real :: gauss_formula_float
   integer :: gauss_formula_int 
   print*, '#This program computes the sum of the first N numbers'
   print*, '#It performs the operation using Integer numbers and FP numbers' 
   print*, '#and prints out the results of both Integer and FP sum'
   print*, '#Operation is performed in two ways:'
   print*, '#  1.summing from 1 to N (from the smallest to the largest)'
   print*, '#  2.summing from N to 1 (from the largest  to the smallest)'  
   print*, '# enter N: '
   read(*,*) nsteps

! direct loop: summing from the smallest one to largest..
! zeroing variables   
   i_rev=0
   directly=0.
   reversly=0.0
   i_dir=0
   do i=1,nsteps
     i_dir = i_dir + i
     ilarge=nsteps +1-i
     i_rev = i_rev + ilarge 
     directly = directly + real(i)
     reversly = reversly + real(ilarge)
   end do 
   k=nsteps
   gauss_formula_float=(real(k)*real(k+1))/2.0 
   gauss_formula_int  =(k*(k+1))/2.0 
    write(*,'(a,i20,g25.15)') '#the mathematical correct result:[n(n+1)/2*]=',gauss_formula_int,gauss_formula_float
    write(*,'(a)') '#summing integers:'
    write(*,'(a,i20,a,i20)') '#direct=', i_dir,'   reverse=',i_rev
    write(*,'(a)') '#summing floating point  numbers' 
    write(*,'(a,i20,a,i20)') '#direct=',int(directly),'   reverse=',int(reversly) 
    write(*,'(a)') '#differences' 
    write(*,'(a,g25.15,a,g25.15)') '#direct=',abs(directly-real(i_dir)),'reverse=',reversly-real(i_rev)


  end program sum_numbers
