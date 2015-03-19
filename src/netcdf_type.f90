module netcdf_type
    implicit none

    integer, parameter :: NO_VALUE = -1000

    type :: netcdf_metadata
      integer     :: ncid
      integer     :: lat_varid, lon_varid, lvl_varid, rec_varid, var_varid
      integer     :: lat_dimid, lon_dimid, lvl_dimid, rec_dimid
      integer     :: lat_n, lon_n, lvl_n, rec_n
      integer     :: ndims
    end type netcdf_metadata

end module netcdf_type
