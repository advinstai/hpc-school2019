interface 
     subroutine greet
     end subroutine
end interface

subroutine greet

  PRINT*,'HELLO, WORLD!'
end subroutine greet

program hello
  call greet
end program
