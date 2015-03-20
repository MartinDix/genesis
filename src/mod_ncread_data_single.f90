module mod_ncread_data_single
  implicit none

contains

subroutine ncread_data_single(nc_meta, var, data, debug)
  use netcdf_type
  use netcdf_check, only: check
  use lower_module, only: to_lower
  use mod_int_to_str, only: int_to_str
  implicit none
  include 'netcdf.inc'

  type(netcdf_metadata), intent(in) :: nc_meta
  character(len=3), intent(in)      :: var
  logical, intent(in)               :: debug
  real, dimension(:), intent(out)   :: data

  integer :: nc_err
  integer :: len_data, len_var
  integer :: varid

  len_data = size(data)

  select case (to_lower(var, 3))
    case ('rec')
      varid = nc_meta%rec_varid
      len_var = nc_meta%rec_n
    case ('lvl')
      varid = nc_meta%lvl_varid
      len_var = nc_meta%lvl_n
    case ('lon')
      varid = nc_meta%lon_varid
      len_var = nc_meta%lon_n
    case ('lat')
      varid = nc_meta%lat_varid
      len_var = nc_meta%lat_n
    case ('var')
      varid = nc_meta%var_varid
      len_var = nc_meta%rec_n * nc_meta%lvl_n * nc_meta%lon_n * nc_meta%lat_n
    case default
      print *, 'expected one of var, lat, lon, lvl, rec, received ' // var
  end select

  if (varid == NO_VALUE) stop "Trying to read a variable that isn't there"
  if (len_data < len_var) then
    print *, "Trying to read " // var // &
    &     " but array isn't big enough."
    print *, "Needed: " // int_to_str(len_var, 8)
    print *, "Given:  " // int_to_str(len_data, 8)
    stop "ERROR"
  end if

  if (debug) print *, "Read " // var

  nc_err = nf_get_var_real(nc_meta%ncid, varid, data)
  call check(nc_err, message = 'in nf_get_var_real for ' // var)

  return
end subroutine ncread_data_single

subroutine ncread_data_single_4D(nc_meta, var, data, debug)
  use netcdf_type
  use netcdf_check, only: check
  use lower_module, only: to_lower
  use mod_int_to_str, only: int_to_str
  implicit none
  include 'netcdf.inc'

  type(netcdf_metadata), intent(in)           :: nc_meta
  character(len=3), intent(in)                :: var
  logical, intent(in)                         :: debug
  real, dimension(:, :, :, :), intent(out)    :: data

  real, dimension(:), allocatable :: tmp_data

  allocate(tmp_data(size(data)))
  ncread_data_single(nc_meta, var, tmp_data, debug)
  data = reshape(tmp_data, shape(data))
  return
end subroutine ncread_data_single_4D


end module mod_ncread_data_single
