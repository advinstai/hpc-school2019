      program main3
c
c  
c
      implicit none
      integer mx
      parameter( mx = 100000 )
      real a(mx)
      real one_over3

      integer i
      ! define 1/3 
      one_over3=1.0/3.0
      ! Initialize the array
      
      do i = 1, mx
         a(i) = float(i)
      end do
      
      ! Average the array
      do i = 1, 10000
      	 call xaver( a, mx ,one_over3)
      end do
      write(*,'(E22.14,2X,E22.14,2X,E22.14)')
     &    a(10000), 20000.0e0,
     &    (a(10000)-20000.0e0)
      end


      subroutine xaver( a, mx,one_over3 )
c
c  Average the data in i direction.
c
      implicit none
      integer mx
      real a(mx)
      real one_over3
      integer i
      
      do i = 1, mx-2
         a(i) = ( a(i) + a(i+1) + a(i+2) ) * one_over3
      end do

      end

