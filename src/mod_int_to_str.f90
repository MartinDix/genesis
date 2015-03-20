module mod_int_to_str
    implicit none
contains
  function int_to_str(i, l) result (s)
    implicit none
    integer, intent(in) :: i
    integer, intent(in) :: l
    character(len=l)    :: s
    character(len=len('(IXX)')) :: f

    write(f, '("(I", I2.2, ")")') l
    write(s, f) i
    return
  end function int_to_str
end module mod_int_to_str
