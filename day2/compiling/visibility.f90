module func
  integer, parameter :: val3=-10, val4=20
  integer :: val5, val6
contains
  integer function add_abs(v1,v2)
    integer, intent(in) :: v1, v2
    add_abs = iabs(v1)+iabs(v2)
  end function add_abs
end module func

program visibility
  use func
  integer :: val1,val2

  val1 = 10
  val2 = -20
  val5 = -5
  val6 = -5
  print*, add_abs(val1,val2), add_abs(val3,val4), add_abs(val5,val6)
end program
  
