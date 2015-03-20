module netcdf_check
  implicit none
  include 'netcdf.inc'

contains
  subroutine check(status, message)
    integer, intent ( in) :: status
    character(len=*), intent(in), optional :: message

    if(status /= nf_noerr) then
      print *, 'A netcdf operation returned an error:'
      print *, trim(nf_strerror(status))
      if (present(message)) then
        print *, message
      end if
      stop "Stopped"
    end if
  end subroutine check
end module netcdf_check
