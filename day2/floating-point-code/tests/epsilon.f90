! compute the machine precision epsilon

program epsilon
  implicit none
  integer, parameter :: prec = 8
  real(kind=prec) :: eps

  eps = 1.0_prec

  print*,'current epsilon, 1 + current epsilon'
  do
    write(6,'(2F20.16)') eps, (1.0_prec + eps)
    eps = eps / 2.0_prec
    if (1.0_prec == 1.0_prec + eps) exit
  end do
  print*,'final epsilon is', eps

end program epsilon
