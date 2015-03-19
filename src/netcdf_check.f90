module netcdf_check
  implicit none
  include 'netcdf.inc'

contains
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf_noerr) then
      print *, trim(nf_strerror(status))
      stop "Stopped"
    end if
  end subroutine check
end module netcdf_check
