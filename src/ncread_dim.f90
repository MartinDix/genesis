subroutine ncread_dim(nc_meta, debug, infile)

  !  First read of NetCDF to get array dimensions
  !  - vjb 23/4/09

  use global, only: num
  use netcdf_check, only: check
  use netcdf_type

  !        use netcdf
  implicit none
  include 'netcdf.inc'

  !------ Passed variables
  !
  logical,intent(in)  :: debug
  character*25,intent(in)  :: infile

  type(netcdf_metadata), intent(out) :: nc_meta

  integer  :: m  ! counters

  !------ NetCDF I/O
  !
  integer  :: ndim
  integer  :: nvar
  
  integer  :: nattr,unlimdimid

  CHARACTER(len=NF_MAX_NAME), DIMENSION(:), ALLOCATABLE :: dim_names
  CHARACTER(len=NF_MAX_NAME), DIMENSION(:), ALLOCATABLE :: var_names
  integer :: i, j
  integer :: var_ndims
  integer, dimension(4) ::var_dimids
  logical :: is_dim

  logical  :: dd

  !------ Open file and fill arrays ------
  !

  dd = .false.

  if (debug) print *,'Reading: ',infile

  call check(nf_open(infile, NF_NOWRITE, nc_meta%ncid))

  call check(nf_inq(nc_meta%ncid,ndim,nvar,nattr,unlimdimid))

  if ((nvar - ndim) <= 0) stop "No Variables found"
  if ((nvar - ndim) > 1) then
    print *, "WARNING, multiple variables possible"
  end if

  print *,'nvar',nvar

  allocate(dim_names(ndim))
  allocate(var_names(nvar))

  dimnames_loop: do i = 1, ndim
    call check(nf_inq_dimname(nc_meta%ncid, i, dim_names(i)))
  end do dimnames_loop

  varnames_loop: do i = 1, nvar
    call check(nf_inq_varname(nc_meta%ncid, i, var_names(i)))
  end do varnames_loop

  nc_meta%var_varid = NO_VALUE
  write(*, 101) 'S', 'Name', '#dims'
  select_var_loop: do i = 1, nvar
    is_dim = .FALSE.
    do j = 1, ndim
      if (var_names(i) == dim_names(j)) is_dim = .TRUE.
    end do
    if (is_dim) then
      write(*, 101) ' ', var_names(i), 'Dim'
    else
      call check(nf_inq_varndims(nc_meta%ncid, i, var_ndims))
      if ((nc_meta%var_varid == NO_VALUE) .and. (var_ndims >= 3) .and. (var_ndims <= 4)) then
        write(*, 100) '*', var_names(i), var_ndims
        nc_meta%var_varid = i
      else
        write(*, 100) ' ', var_names(i), var_ndims
      end if
    end if ! is_dim
  end do select_var_loop
100 FORMAT(1X, A1, 1X, A25, 1X, I4) ! For vars
101 FORMAT(1X, A1, 1X, A25, 1X, A4) ! For header and dims

  if (nc_meta%var_varid < 0) stop "ERROR: NO VARIABLE FOUND"

  call check(nf_inq_varndims(nc_meta%ncid, nc_meta%var_varid, nc_meta%ndims))
  call check(nf_inq_vardimid(nc_meta%ncid, nc_meta%var_varid, var_dimids))

  dd = (nc_meta%ndims == 3)
  if (dd) then
    var_dimids(4) = var_dimids(3)
    var_dimids(3) = -1000
  end if
  nc_meta%rec_dimid = var_dimids(4)
  nc_meta%lvl_dimid = var_dimids(3)
  nc_meta%lat_dimid = var_dimids(2)
  nc_meta%lon_dimid = var_dimids(1)
  do j = 1, nvar
    if (var_names(j) == dim_names(nc_meta%rec_dimid)) nc_meta%rec_varid = j
    if (.not. dd) if (var_names(j) == dim_names(nc_meta%lvl_dimid)) nc_meta%lvl_varid = j
    if (var_names(j) == dim_names(nc_meta%lon_dimid)) nc_meta%lon_varid = j
    if (var_names(j) == dim_names(nc_meta%lat_dimid)) nc_meta%lat_varid = j
  end do

  print *, 'assuming the following dimensions:'
  print *, 'lon: ' // trim(var_names(nc_meta%lon_varid))
  print *, 'lat: ' // trim(var_names(nc_meta%lat_varid))
  if (.not. dd) print *, 'lvl: ' // trim(var_names(nc_meta%lvl_varid))
  print *, 'rec: ' // trim(var_names(nc_meta%rec_varid))

  if (debug) print*,'..1..'

  call check(nf_inq_dimlen(nc_meta%ncid, nc_meta%lon_varid, nc_meta%lon_n))
  call check(nf_inq_dimlen(nc_meta%ncid, nc_meta%lat_varid, nc_meta%lat_n))
  if (dd)then ! Reduced dimensions
    nc_meta%lvl_n = 1
  else
    call check(nf_inq_dimlen(nc_meta%ncid, nc_meta%lvl_varid, nc_meta%lvl_n))
  endif
  call check(nf_inq_dimlen(nc_meta%ncid, nc_meta%rec_varid, nc_meta%rec_n))

  if (debug) print*,'..2..'

  print *,'nlon = ',nc_meta%lon_n
  print *,'nlat = ',nc_meta%lat_n
  print *,'nlvl = ',nc_meta%lvl_n
  print *,'nrec = ',nc_meta%rec_n

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
  !-------------------------------------------
  !  flip 1d array in vertical

  !  subroutine arrswapk(nx,ny,nz,nc_meta%rec_n,infld)
  subroutine arrswapk(nz,infld)

    integer nz
    !  real, intent(inout)  :: infld(nx,ny,nz,nc_meta%rec_n)
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
  !  flip 3d array in vertical - repeat for nc_meta%rec_ns

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
