module global

  implicit none

  save

  integer,parameter  :: num = 20
  integer,parameter  :: max_nrecs = 500
  integer,parameter  :: NDIMS = 4
  integer,parameter  :: NMSL  = 1
  integer,parameter  :: nqlev = 12
  integer,parameter  :: maxfiles = 7
  integer,parameter  :: nvars = 7
  integer,parameter  :: acc = 0.007  ! lat/lon acc
  real,parameter  :: gravity = 9.80665
  real,parameter  :: rho = 1.225
  real,parameter  :: minp = 0.0
  real,parameter  :: maxz = 43000.0
  real,parameter  :: secday = 86400.0
  real,parameter  :: omg = 7.292e-05
  real,parameter  :: rcp = 287./1004.
  real,parameter  :: pi = 3.141592654     ! 4*atan(1.0)
  real,parameter  :: zero = 0.0
  real,parameter  :: zrough = 0.02  ! z0 grass
  real,parameter  :: vkman = 0.4  ! Von Karman const
  integer,parameter  :: ntmax=50000
  integer,parameter  :: nzmax=100
  integer,parameter  :: umlev=38  ! UM 38 levels

end module global
