	subroutine geostr(frclat,frclon,nrec,nlev,nlat,nlon,levs,lats,lons,  &
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
	real,dimension(nrec),intent(in)		:: frclat,frclon
	integer,intent(in)			:: nrec,nlev,nlat,nlon
	real,dimension(nlev,nrec),intent(in)	:: t
	
	real,dimension(nlon,nlat,nrec,nlev),intent(in)	:: z_in
	
	real,dimension(nlev),intent(in)		:: levs
	real,dimension(nlat),intent(in)		:: lats
	real,dimension(nlon),intent(in)		:: lons
	integer,dimension(nrec),intent(in)	:: fynindex,fysindex
	integer,dimension(nrec),intent(in)	:: fxeindex,fxwindex
	real,dimension(nrec),intent(in)		:: dx,dy
	logical,intent(in)			:: debug

!------ Working ------
!
	integer		                        :: ynindex,ysindex
	integer		                        :: xeindex,xwindex
	integer					:: k,rec,i
!	real,dimension(:,:),allocatable		:: pt
	real,dimension(:,:),allocatable		:: dphidx,dphidy
	real					:: radlat
	
	real					:: surfdist
	real					:: radian
	
!------ Output ------
!
	real,dimension(nlev,nrec),intent(out)		:: ug,vg	
	real,dimension(nlev,nrec),intent(out)		:: pt
	real,intent(out)				:: f

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
	 radlat = radian(frclat(rec))		! use the radian function
	 f = 2*omg*sin(radlat)			! Coriolis term


	 do k=1,nlev
     
	  dphidx(k,rec) = (z_in(fxeindex(rec),fynindex(rec),k,rec)              &
     &                     -z_in(fxwindex(rec),fynindex(rec),k,rec))/dx(rec)
     
	  dphidy(k,rec) = (z_in(fxeindex(rec),fynindex(rec),k,rec)              &
     &                     -z_in(fxeindex(rec),fysindex(rec),k,rec))/dy(rec)
	  
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
