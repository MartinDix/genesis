subroutine find_latlon(nrecs,nlat,nlon,finlat,finlon,lat,lon,fynindex,fysindex, &
    &				fxeindex,fxwindex,debug)

    !	From user defined lat/lon - finds nearest neighbour lat/lon on
    !	grid array (for linear interpolation)
    !	- vjb 2/2/09

    use global

    implicit none

    real		                        :: inlat,inlon
    integer		                        :: ynindex,ysindex
    integer		                        :: xeindex,xwindex
    real,dimension(nrecs),intent(in)		:: finlat,finlon

    integer,intent(in)			:: nlat,nlon,nrecs

    real,dimension(nlat),intent(in)		:: lat
    real,dimension(nlon),intent(in)		:: lon

    integer,dimension(nrecs),intent(out)	:: fynindex,fysindex
    integer,dimension(nrecs),intent(out)	:: fxeindex,fxwindex

    integer					:: i
    logical					:: same_y,same_x

    logical					:: positive	! search up/down
    logical,intent(in)			:: debug


    do i=1,nrecs

        inlat=finlat(i)
        inlon=finlon(i)
        same_x = .false.
        same_y = .false.

        !	if (debug) print *,'nlat,nlon,inlat,inlon = ',nlat,nlon,inlat,inlon


        !---------------- Find nearest lat/lon interval -------------
        !  First, the nearest northern latitude

        positive = .true.
        !	print *,'inlat=',inlat
        !	print *,'lat=',lat
        !	print *,'nlat=',nlat
        call search_grid(inlat,lat,nlat,positive,ynindex,debug)
        !	if (debug) print *,'ynindex,inlat,lat(ynindex) = ',ynindex,inlat,lat(ynindex)


        !------------------------------------------------------------
        !  Now, do the same in southerly direction to get nearest south latitude

        positive = .false.
        call search_grid(inlat,lat,nlat,positive,ysindex,debug)
        !	if (debug) print *,'ysindex,inlat,lat(ysindex) = ',ysindex,inlat,lat(ysindex)


        !------------------------------------------------------------
        !  Once again, now scanning west to east...

        positive = .true.
        call search_grid(inlon,lon,nlon,positive,xeindex,debug)
        !	if (debug) print *,'xeindex,inlon,lon(xeindex) = ',xeindex,inlon,lon(xeindex)


        !----------------------------------------------------------
        !  Finally, scan east to west...

        positive = .false.
        call search_grid(inlon,lon,nlon,positive,xwindex,debug)
        !	if (debug) print *,'xwindex,inlon,lon(xwindex) = ',xwindex,inlon,lon(xwindex)


        !	print *,'Found nearest gridpoints to specified lat/lon'
        print *,'glrindex',ynindex,ysindex,xeindex,xwindex
        print *,'glrindex',lat(ynindex),lat(ysindex),lon(xeindex),lon(xwindex)

        !--------- Check to see if indices are the same --------
        !
        if (ynindex.eq.ysindex)then
            same_y = .true.
        endif
        if (xeindex.eq.xwindex)then
            same_x = .true.
        endif

        if (same_y.and.debug) print *,'N and S latitudes are the same: Y'
        if (same_x.and.debug) print *,'E and W latitudes are the same: X'

        fynindex(i)=ynindex
        fysindex(i)=ysindex
        fxeindex(i)=xeindex
        fxwindex(i)=xwindex

    enddo

    return

end subroutine find_latlon
