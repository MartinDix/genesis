subroutine ncread_dim(ncid,nvar,nlon,nlat,nlev,nrec,  &
  &                       debug,infile, varid, var_dimids)

  !  First read of NetCDF to get array dimensions
  !  - vjb 23/4/09

  use global, only: num

  !        use netcdf
  implicit none
  include 'netcdf.inc'

  !------ Passed variables
  !
  logical,intent(in)  :: debug
  character*25,intent(in)  :: infile

  integer,intent(out)  :: nlon  !
  integer,intent(out)  :: nlat  ! passed
  integer,intent(out)  :: nlev  ! arrays 2
  integer,intent(out)  :: nrec  !

  integer  :: m  ! counters

  !------ NetCDF I/O
  !
  integer,intent(out)  :: ncid
  integer  :: ndim
  integer,intent(out)  :: nvar
  
  integer  :: nattr,unlimdimid

  integer  :: lat_varid
  integer  :: lon_varid
  integer  :: lvl_varid
  integer  :: rec_varid
  

  CHARACTER(len=NF_MAX_NAME), DIMENSION(:), ALLOCATABLE :: dim_names
  CHARACTER(len=NF_MAX_NAME), DIMENSION(:), ALLOCATABLE :: var_names
  integer :: i, j
  integer, intent(out) :: varid
  integer :: var_ndims
  integer, dimension(4), intent(out) ::var_dimids
  logical :: is_dim

  logical  :: dd

  !------ Open file and fill arrays ------
  !

  dd = .false.

  if (debug) print *,'Reading: ',infile

  call check(nf_open(infile,NF_NOWRITE,ncid))

  call check(nf_inq(ncid,ndim,nvar,nattr,unlimdimid))

  if ((nvar - ndim) <= 0) stop "No Variables found"
  if ((nvar - ndim) > 1) then
    print *, "WARNING, multiple variables possible"
  end if

  print *,'nvar',nvar

  allocate(dim_names(ndim))
  allocate(var_names(nvar))

  dimnames_loop: do i = 1, ndim
    call check(nf_inq_dimname(ncid, i, dim_names(i)))
  end do dimnames_loop

  varnames_loop: do i = 1, nvar
    call check(nf_inq_varname(ncid, i, var_names(i)))
  end do varnames_loop

  varid = -1000
  write(*, 101) 'S', 'Name', '#dims'
  select_var_loop: do i = 1, nvar
    is_dim = .FALSE.
    do j = 1, ndim
      if (var_names(i) == dim_names(j)) is_dim = .TRUE.
    end do
    if (is_dim) then
      write(*, 101) ' ', var_names(i), 'Dim'
    else
      call check(nf_inq_varndims(ncid, i, var_ndims))
      if ((varid < 0) .and. (var_ndims >= 3) .and. (var_ndims <= 4)) then
        write(*, 100) '*', var_names(i), var_ndims
        varid = i
      else
        write(*, 100) ' ', var_names(i), var_ndims
      end if
    end if ! is_dim
  end do select_var_loop
100 FORMAT(1X, A1, 1X, A25, 1X, I4) ! For vars
101 FORMAT(1X, A1, 1X, A25, 1X, A4) ! For header and dims

  if (varid < 0) stop "ERROR: NO VARIABLE FOUND"

  call check(nf_inq_varndims(ncid, varid, var_ndims))
  call check(nf_inq_vardimid(ncid, varid, var_dimids))

  dd = (var_ndims == 3)
  if (dd) then
    var_dimids(4) = var_dimids(3)
    var_dimids(3) = -1000
  end if
  rec_varid = var_dimids(1)
  lvl_varid = var_dimids(2)
  lat_varid = var_dimids(3)
  lon_varid = var_dimids(4)

  print *, 'assuming the following dimensions:'
  print *, 'lon: ' // trim(dim_names(lon_varid))
  print *, 'lat: ' // trim(dim_names(lat_varid))
  if (.not. dd) print *, 'lvl: ' // trim(dim_names(lvl_varid))
  print *, 'rec: ' // trim(dim_names(rec_varid))

  if (debug) print*,'..1..'

  call check(nf_inq_dimlen(ncid, lon_varid, nlon))
  call check(nf_inq_dimlen(ncid, lat_varid, nlat))
  if (dd)then ! Reduced dimensions
    nlev = 1
  else
    call check(nf_inq_dimlen(ncid, lvl_varid, nlev))
  endif
  call check(nf_inq_dimlen(ncid, rec_varid, nrec))

  if (debug) print*,'..2..'

  print *,'nlon = ',nlon
  print *,'nlat = ',nlat
  print *,'nlev = ',nlev
  print *,'nrec = ',nrec

  dd = .false.

  !------- return
  !
  print *,"**** Successfully read input file: ",infile

  ! For compatibility:
  nvar = var_ndims + 1

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

  !-------------------------------------------
  !  flip 1d array in vertical

  !  subroutine arrswapk(nx,ny,nz,nrec,infld)
  subroutine arrswapk(nz,infld)

    integer nz
    !  real, intent(inout)  :: infld(nx,ny,nz,nrec)
    real, intent(inout)  :: infld(nz)
    real  :: outfld(nz)

    integer k
    integer kinv
      


    do k=1,nz
      kinv = nz+1 - k
      outfld(k) = infld(kinv)
    enddo
    do k=1,nz
      infld(k) = outfld(k)
    enddo

    return

  end subroutine arrswapk

  !----------------------------------------------
  !  flip 3d array in vertical - repeat for nrecs

  subroutine arrswap3dk(nx,ny,nz,nrec,infld)

    integer   :: nx,ny,nz,nrec
    real, intent(inout)  :: infld(nx,ny,nz,nrec)
    real  :: outfld(nx,ny,nz,nrec)

    integer  :: j,k,l
    integer  :: kinv
      

    do j=1,nrec
      do k=1,nz
        kinv = nz+1 - k
        do l=1,ny
          do m=1,nx
            outfld(m,l,k,j) = infld(m,l,kinv,j)
          enddo
        enddo
      enddo
    enddo

    do j=1,nrec
      do k=1,nz
        do l=1,ny
          do m=1,nx
            infld(m,l,k,j) = outfld(m,l,k,j)
          enddo
        enddo
      enddo
    enddo

    return

  end subroutine arrswap3dk

  !------------------------------------------

end subroutine ncread_dim
