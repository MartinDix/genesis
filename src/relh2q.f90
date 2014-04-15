
	subroutine relh2q(nlon,nlat,nlev,nrec,lev,temp,hum)

!	Takes input temperature (K) and relative humidity (%) and calculates corresponding
!	specific humidity (kg/kg)
!
!	- vjb 31/1/2011

	implicit none

!------ Declarations ------
!
	integer,intent(in)					:: nlon
	integer,intent(in)					:: nlat
	integer,intent(in)					:: nlev
	integer,intent(in)					:: nrec

	real,dimension(nlev),intent(in)				:: lev
	real,dimension(nlon,nlat,nlev,nrec),intent(in)		:: temp   ! temperature
	real							:: temp1  ! temp T

	real,dimension(nlon,nlat,nlev,nrec),intent(inout)	:: hum    ! relh => q
	real							:: es1    ! sat. vap. press.
	real							:: vap    ! vap. press.

	integer							:: i,j,k,m


!------ Kick off ------
!

	do m=1,nrec
	  do k=1,nlev
	    do j=1,nlat
	      do i=1,nlon

                temp1 = temp(i,j,k,m)
	        call establ(es1,temp1)
                vap = es1*(hum(i,j,k,m)/100)
                hum(i,j,k,m) = (vap*0.622)/(lev(k)*100)

	      enddo
	    enddo
	  enddo
	enddo

	return

	end subroutine relh2q
