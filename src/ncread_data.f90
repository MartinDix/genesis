subroutine ncread_data(nc_meta, lon,lat,lev,rec,dat)

  !  fills arrays according to NetCDF dimensions
  !  - vjb 20/8/09

  use netcdf_check, only: check
  use netcdf_type

  !        use netcdf
  implicit none
  include 'netcdf.inc'
  
  type(netcdf_metadata), intent(in)  :: nc_meta

  real,dimension(nc_meta%lon_n),intent(inout)  :: lon
  real,dimension(nc_meta%lat_n),intent(inout)  :: lat
  real,dimension(nc_meta%lvl_n),intent(inout)  :: lev
  real,dimension(nc_meta%rec_n),intent(inout)  :: rec

  real,dimension(nc_meta%lon_n,nc_meta%lat_n,nc_meta%lvl_n,nc_meta%rec_n),intent(inout)  :: dat

  integer  :: lat_varid
  integer  :: lon_varid
  integer  :: lvl_varid
  integer  :: rec_varid
  integer  :: dat_varid

  logical  :: dd

  print *,'ncread_data: shape(dat) = ',shape(dat)

  rec_varid = nc_meta%rec_varid
  lvl_varid = nc_meta%lvl_varid
  lat_varid = nc_meta%lat_varid
  lon_varid = nc_meta%lon_varid
  dat_varid = nc_meta%var_varid

  dd = (nc_meta%ndims == 3)

  call check(nf_get_var_real(nc_meta%ncid,lon_varid,lon))
  call check(nf_get_var_real(nc_meta%ncid,lat_varid,lat))
  call check(nf_get_var_real(nc_meta%ncid,rec_varid,rec))
  if (.not. dd) then
    call check(nf_get_var_real(nc_meta%ncid,lvl_varid,lev))
  else
    lev = 1000.0
  end if

  call check(nf_get_var_real(nc_meta%ncid,dat_varid,dat))

  call check(nf_close(nc_meta%ncid))

  return

  !------ Subroutines ------
  !  NetCDF function

end subroutine ncread_data
