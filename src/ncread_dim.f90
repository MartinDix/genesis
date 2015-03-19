subroutine ncread_dim(ncid,nvar,nlon,nlat,nlev,nrec,  &
  &                       debug,infile)

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
  integer  :: dat_varid

  logical  :: dd

  !------ Open file and fill arrays ------
  !

  dd = .false.

  if (debug) print *,'Reading: ',infile

  call check(nf_open(infile,NF_NOWRITE,ncid))

  call check(nf_inq(ncid,ndim,nvar,nattr,unlimdimid))

  print *,'nvar',nvar

  lon_varid = 1
  lat_varid = 2
  if (nvar.eq.4)then
    rec_varid = 3  ! problem with ncid
    dat_varid = 4
  elseif(nvar.eq.5)then
    lvl_varid = 3
    rec_varid = 4
    dat_varid = 5
  else
    stop 'ERROR: nvar dimension outside acceptable threshold'
  endif

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
