program loop
  implicit none
  integer :: i
  real    :: x

  i = 0
  x = 0.0
  do 
    i = i + 1
    x = x + 0.01
    if (x >= 1.0) exit
  end do
  print*, 'Number of iterations using real counter:   ', i

  i = 0
  x = 0.0
  do
    i = i + 1
    x = x + 0.01
    if (i >= 100) exit
  end do
  print*, 'Number of iterations using integer counter:', i

  print*, 'Difference 100*0.01 vs. Sum_100(0.01): ', 1.0-x

end program loop
