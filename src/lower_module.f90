module lower_module
  implicit none

contains
  function to_lower(string, string_len) result (lower_string)
    implicit none
    integer, intent(in) :: string_len
    character(len=string_len), intent(in) :: string
    character(len=string_len) :: lower_string

    integer :: i, ic

    lower_string = string

    do i = 1, string_len
      ic = ichar(string(i:i))
      if (ic >= ichar('A') .and. ic <= ichar('Z')) then
        lower_string(i:i) = char(ic+32)
      end if
    end do
    return
  end function to_lower

end module lower_module
