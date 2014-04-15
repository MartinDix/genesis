
	subroutine search_grid(usrpt,dat,max,positive,index,debug)

!	Given a grid point (usrpt), will find nearest neighbour in dat array.
!	Choose positive = .true. to search 'up' and .false. to search 'down'.

	implicit none

!------- declarations 

	real,parameter			:: acc = 0.007 ! accuracy of real number search
	integer				:: max ! size of destination array
	real,intent(in)			:: usrpt ! user defined grid point (f5.1)
	real,dimension(max),intent(in)	:: dat ! destination grid array
	logical				:: positive ! search up or down (+/-)
	integer,intent(out)		:: index ! output grid index
	logical				:: debug ! debug option
        real				:: modulus,dd,delta,degsearch,gacc
	integer				:: mmax,nn

!------- working

	integer				:: m,n	! counters
	real				:: gridpt ! working grid point value

!---------------- Search -------------
!

	degsearch = 5.0 !valid for access-g
	degsearch = 0.5 !valid for access-c where dlat=0.05deg =>steps of 0.001
        mmax = 500
        delta = degsearch/float(mmax)
        gacc = 0.5*delta
	if (positive)then

	 gridpt = usrpt				! initialise
	 do m=1,mmax				! search 'up' 5 deg in steps of 0.01
	  do n=1,max
	   dd = gridpt-dat(n)
	   if (modulus(dd).le.gacc)then
	    goto 11
	   endif
	  enddo
	  gridpt = gridpt+delta
	 enddo

11	 continue
	
	 index = n

	else					! search 'down' 5 deg in steps of 0.01

	 gridpt = usrpt
	 do m=1,mmax
	  do n=max,1,-1
	   dd = gridpt-dat(n)
	   if (modulus(dd).le.gacc)then
	    goto 12
	   endif
	  enddo
	  gridpt = gridpt-delta
	 enddo

12	 continue

	 index = n

	endif
	
!------------------------------------------------------------
!my search calc
!	 gridpt = usrpt				! initialise
!          nn=-1
!	  do n=1,max
!	   dd = gridpt-dat(n)
!           if(nn.lt.0.and.dd.le.0.0) nn = n
!           print *,'n,nn,dd,dat(n),usrpt=',n,nn,dd,dat(n),usrpt
!	  enddo
!	if (positive)then
!	  index = nn
!	else
!	  index = nn - 1
!	endif
!------------------------------------------------------------

	return

	end
