program testsum
  implicit none

  integer, parameter :: n=200
  integer :: data(n), asum, i

  do i=1,200
    data(i) = i-100
  end do

  call sum_abs(data,n,asum)
  print*, 'sum=',asum
end program testsum 
