subroutine geostr(frclat,nrec,nlev,nlat,nlon,levs,  &
  &                    fynindex,fysindex,fxeindex,fxwindex,t,z_in,dx,dy,ug,vg,f,pt,   &
  &                    debug)

  !------- Derive geostrophic wind profiles --------
  !
  !	Find Coriolis term and use the grid geopotential
  !	finite differences to get d(phi)/dx and d(phi)/dy
  !	(note: using pressure coordinates).
  !

  use global

  implicit none

  !------ Input ------
  !
  real,dimension(nrec),intent(in)       :: frclat
  integer,intent(in)                    :: nrec,nlev,nlat,nlon
  real,dimension(nlev,nrec),intent(in)  :: t

  !	real,dimension(nlon,nlat,nrec,nlev),intent(in)	:: z_in
  real,dimension(nlon,nlat,nlev,nrec),intent(in)    :: z_in

  real,dimension(nlev),intent(in)     :: levs
  integer,dimension(nrec),intent(in)  :: fynindex,fysindex
  integer,dimension(nrec),intent(in)  :: fxeindex,fxwindex
  real,dimension(nrec),intent(in)     :: dx,dy
  logical,intent(in)                  :: debug

  !------ Working ------
  !
  real :: radian  ! function declared in radian.f90
  
  integer   :: k,rec
  !	real,dimension(:,:),allocatable		:: pt
  real,dimension(:,:),allocatable   :: dphidx,dphidy
  real          :: radlat

  
  !------ Output ------
  !
  real,dimension(nlev,nrec),intent(out)   :: ug,vg
  real,dimension(nlev,nrec),intent(out)   :: pt
  real,intent(out)                        :: f

  !------- Create potential temperature profiles -------
  !
  print *,'Deriving potential temperatures'

  do rec=1,nrec
    !	do rec=2,2
    do k=1,nlev
      pt(k,rec) = t(k,rec)*((100000/levs(k))**rcp)
    enddo
  enddo

  !	print *,'t = ',t
  !	if (debug) print *,'t,pt = ',(t(k,2),k=1,nlev),(pt(k,2),k=1,nlev)

  !------- Derive geostrophic wind profiles --------
  !
  !	Find Coriolis term and use the grid geopotential
  !	finite differences to get d(phi)/dx and d(phi)/dy
  !	(in pressure coordinates).

  print *,'Deriving geostrophic winds'

  allocate(dphidx(nlev,nrec))
  allocate(dphidy(nlev,nrec))

  do rec=1,nrec
    !     	do rec=2,2
    radlat = radian(frclat(rec))  ! use the radian function
    f = 2*omg*sin(radlat)         ! Coriolis term


    do k=1,nlev

      if (dx(rec).eq.0.)then
        dphidx(k,rec) = 0.0
      else
        dphidx(k,rec) = (z_in(fxeindex(rec),fynindex(rec),k,rec)           &
          &                       -z_in(fxwindex(rec),fynindex(rec),k,rec))/dx(rec)
      endif
      if (dy(rec).eq.0.)then
        dphidy(k,rec) = 0.0
      else
        dphidy(k,rec) = (z_in(fxeindex(rec),fynindex(rec),k,rec)           &
          &                       -z_in(fxeindex(rec),fysindex(rec),k,rec))/dy(rec)
      endif

      ug(k,rec) = (-1/f)*(dphidy(k,rec))

      vg(k,rec) = (1/f)*(dphidx(k,rec))

      !	  if (debug) print *,'dx,dy = ',dx(rec),dy(rec)
      !	  if (debug) print *,'ug,vg = ',ug(k,1),vg(k,1)

    enddo

  enddo


  !	if (debug) print *,'dx,dy = ',dx(rec),dy(rec)

  if (debug) print *,'done!'

  return

end subroutine geostr
