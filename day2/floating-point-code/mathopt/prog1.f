      program main1
c
c
      implicit none
      integer mx
      parameter( mx = 100000 )
      real*4 a(mx)

      integer i

      ! Initialize the array
      
      do i = 1, mx
        a(i) = float(i)
      end do
      
      ! Average the array
      do i = 1, 10000
        call xaver( a, mx )
      end do
      write(*,'(E22.14,2X,E22.14,2X,E22.14)')
     &    a(10000), 15000.0e0,
     &    (a(10000)-15000.0e0)
      end


      subroutine xaver( a, mx )
c
c  Average the data in i direction.
c
      implicit none
      integer mx
      real*4 a(mx)
      integer i
      
      do i = 1, mx-1
         a(i) = ( a(i) + a(i+1) ) / 2.0e0
      end do

      end

