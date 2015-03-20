module netcdf_type
    implicit none

    integer, parameter :: NO_VALUE = -1000

    type :: netcdf_metadata
      integer     :: ncid = NO_VALUE
      integer     :: lat_varid = NO_VALUE
      integer     :: lon_varid = NO_VALUE
      integer     :: lvl_varid = NO_VALUE
      integer     :: rec_varid = NO_VALUE
      integer     :: var_varid = NO_VALUE
      integer     :: lat_dimid = NO_VALUE
      integer     :: lon_dimid = NO_VALUE
      integer     :: lvl_dimid = NO_VALUE
      integer     :: rec_dimid = NO_VALUE
      integer     :: lat_n = NO_VALUE
      integer     :: lon_n = NO_VALUE
      integer     :: lvl_n = NO_VALUE
      integer     :: rec_n = NO_VALUE
      integer     :: ndims = NO_VALUE
    end type netcdf_metadata

end module netcdf_type
