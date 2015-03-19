subroutine ncread_data(ncid,nvar,nlon,nlat,nlev,nrec,  &
  &                          lon,lat,lev,rec,dat, varid, var_dimids)

  !  fills arrays according to NetCDF dimensions
  !  - vjb 20/8/09

  use global

  !        use netcdf
  implicit none
  include 'netcdf.inc'

  integer  :: i
  
  integer,intent(in)  :: ncid
  integer,intent(in)  :: nvar

  integer,intent(in)  :: nlon
  integer,intent(in)  :: nlat
  integer,intent(in)  :: nlev
  integer,intent(in)  :: nrec

  integer, intent(in) :: varid
  integer, dimension(4), intent(in) :: var_dimids

  real,dimension(nlon),intent(inout)  :: lon
  real,dimension(nlat),intent(inout)  :: lat
  real,dimension(nlev),intent(inout)  :: lev
  real,dimension(nrec),intent(inout)  :: rec

  real,dimension(nlon,nlat,nlev,nrec),intent(inout)  :: dat

  character*20  :: var(num)
  character*20  :: name

  integer  :: lat_varid
  integer  :: lon_varid
  integer  :: lvl_varid
  integer  :: rec_varid
  integer  :: dat_varid

  logical                                         :: dd


  dd = .false.

  print *,'ncread_data: shape(dat) = ',shape(dat)

  rec_varid = var_dimids(1)
  lvl_varid = var_dimids(2)
  lat_varid = var_dimids(3)
  lon_varid = var_dimids(4)
  dat_varid = varid

  dd = (lvl_varid < 0)

  write(*, *) var_dimids, varid, dd

  call check(nf_get_var_real(ncid,lon_varid,lon))
  call check(nf_get_var_real(ncid,lat_varid,lat))
  call check(nf_get_var_real(ncid,rec_varid,rec))
  if (.not. dd) then
    call check(nf_get_var_real(ncid,lvl_varid,lev))
  else
    lev = 1000.0
  end if

  call check(nf_get_var_real(ncid,dat_varid,dat))

  call check(nf_close(ncid))

  return

  !------ Subroutines ------
  !  NetCDF function

contains

  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf_noerr) then
      print *, trim(nf_strerror(status))
      stop "Stopped"
    end if
  end subroutine check

end subroutine ncread_data
