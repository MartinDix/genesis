subroutine interp_linear(nlat,nlon,nlev,nrec,fynindex,fysindex, &
    &                           fxeindex,fxwindex,frclat,frclon,lats,  &
    &                           lons,levs,same_x,same_y,dat1,dat2,dat3,   &
    &                           dat4,dat5,dat6,dat7,i_int1,i_int2,i_int3, &
    &                           i_int4,i_int5,i_int6,i_int7,dx,dy,gradt,gradq,    &
    &                           debug)

    !	Takes nearest neighbour lat/lon to user defined lat/lon and linearly
    !	interpolates input fields at all levels
    !	- vjb 2/2/2009
    !
    !       Swapped array dimensions around for all *_in arrays to match input dimensions
    !       - vjb 10/11/2010
    !
    !	Re-check of array dimensions with addition of ncread_data to main program
    !	and space set out for incorporation of vertical advective tendencies
    !	- vjb 16/12/2010

    use global

    implicit none

    !------ Input -----
    !
    integer,intent(in)				:: nlat,nlon,nlev,nrec
    integer,dimension(nrec),intent(in)		:: fynindex,fysindex
    integer,dimension(nrec),intent(in)		:: fxeindex,fxwindex

    real,dimension(nrec),intent(in)			:: frclat,frclon
    real,dimension(nlat),intent(in)			:: lats
    real,dimension(nlon),intent(in)			:: lons
    real,dimension(nlev),intent(in)			:: levs

    real,dimension(nlon,nlat,NMSL,nrec),intent(in)	:: dat1
    real,dimension(nlon,nlat,nlev,nrec),intent(in)	:: dat2,dat3
    real,dimension(nlon,nlat,nlev,nrec),intent(in)	:: dat4,dat5,dat6,dat7     ! vjb 2.3

    logical,dimension(nrec),intent(in)		:: same_x,same_y
    logical,intent(in)				:: debug

    !------ Working variables/arrays ------
    !
    integer						:: i,j,k,l,m,rec

    real                 				:: surfdist
    real						:: factx,facty

    real,dimension(:,:,:,:),allocatable		:: var_in

    real,dimension(:),allocatable			:: z1_2d,z2_2d
    real,dimension(:),allocatable			:: n_2d,s_2d
    real,dimension(:),allocatable			:: e_2d,w_2d

    real,dimension(:,:),allocatable			:: z1_var,z2_var,i_var
    real,dimension(:,:),allocatable			:: n_var,s_var
    real,dimension(:,:),allocatable			:: e_var,w_var
    real,dimension(nlev,nrec)			:: gradtx,gradty,gradtz
    real,dimension(nlev,nrec)			:: gradqx,gradqy,gradqz

    !------ Output arrays ------
    !
    real,dimension(nrec),intent(out)		:: i_int1
    real,dimension(nlev,nrec),intent(out)		:: i_int2,i_int3,i_int4
    real,dimension(nlev,nrec),intent(out)		:: i_int5,i_int6,i_int7      ! vjb 2.3

    real,dimension(nlev,nrec),intent(out)		:: gradt,gradq
    real,dimension(nrec),intent(out)		:: dx,dy

    !------ Prepare working arrays ------
    !
    allocate(n_2d(nrec))			! north
    allocate(s_2d(nrec))			! south
    allocate(e_2d(nrec))			! east
    allocate(w_2d(nrec))			! west
    allocate(z1_2d(nrec))			! interpolated N/S points
    allocate(z2_2d(nrec))			! interpolated E/W points

    allocate(z1_var(nlev,nrec))		! interpolated N/S points
    allocate(z2_var(nlev,nrec))		! interpolated E/W points
    allocate(n_var(nlev,nrec))		! north
    allocate(s_var(nlev,nrec))		! south
    allocate(e_var(nlev,nrec))		! east
    allocate(w_var(nlev,nrec))		! west
    allocate(i_var(nlev,nrec))		! final interpolated output

    allocate(var_in(nlon,nlat,nlev,nrec))	! working array


    !	print *,'dat1 = ',dat1(1,1,1,2)
    !	print *,'dat2 = ',dat2(64,65,10,3)
    !	print *,'dat3 = ',dat3(1,1,1,4)
    !	print *,'dat4 = ',dat4(1,1,1,1)
    !	print *,'dat5 = ',dat5(1,1,1,1)
    !	print *,'dat6 = ',dat6(20,20,10,1)
    !        print *,'dat7 = ',dat7(20,20,10,1)                            ! vjb 2.3



    !------ Loop over input lonlat values -----
    !
    do rec=1,nrec
        !	do rec=1,1

        !------ Calculate dx and dy for later gradient calculations -----
        !
        !         print *,'glr11aA shape2 nlev,nrec',nlev,nrec
        !         print *,'glr11aA shape2',shape(gradt),shape(gradq)

        if((fxeindex(rec).eq.fxwindex(rec)).and.(fynindex(rec).eq.fysindex(rec)))then

            dx(rec) = surfdist(lats(fynindex(rec)),lons(fxeindex(rec)+1),lats(fynindex(rec)),    &
                &	        lons(fxwindex(rec)))
            dy(rec) = surfdist(lats(fynindex(rec)+1),lons(fxeindex(rec)),lats(fysindex(rec)),    &
                &	        lons(fxeindex(rec)))

        elseif (fxeindex(rec).eq.fxwindex(rec))then

            dx(rec) = surfdist(lats(fynindex(rec)),lons(fxeindex(rec)+1),lats(fynindex(rec)),    &
                &	        lons(fxwindex(rec)))
            dy(rec) = surfdist(lats(fynindex(rec)),lons(fxeindex(rec)),lats(fysindex(rec)),    &
                &	        lons(fxeindex(rec)))

        elseif(fynindex(rec).eq.fysindex(rec))then

            dx(rec) = surfdist(lats(fynindex(rec)),lons(fxeindex(rec)),lats(fynindex(rec)),    &
                &	        lons(fxwindex(rec)))
            dy(rec) = surfdist(lats(fynindex(rec)+1),lons(fxeindex(rec)),lats(fysindex(rec)),    &
                &	        lons(fxeindex(rec)))

        else

            dx(rec) = surfdist(lats(fynindex(rec)),lons(fxeindex(rec)),lats(fynindex(rec)),    &
                &	        lons(fxwindex(rec)))
            dy(rec) = surfdist(lats(fynindex(rec)),lons(fxeindex(rec)),lats(fysindex(rec)),    &
                &	        lons(fxeindex(rec)))

        endif

        !         print *,'glr interp.. rec,dx(rec),dy(rec)=',rec,dx(rec),dy(rec)

        if (debug) print *,'frclat(rec),lons(e),lons(w),dx = ', &
            &           frclat(rec),lons(fxeindex(rec)),lons(fxwindex(rec)),dx(rec)

        if (debug) print *,fxeindex(rec),fxwindex(rec)
        if (debug) print *,fynindex(rec),fysindex(rec)

        if (debug) print *,'frclon(rec),lats(n),lats(s),dy = ', &
            &           frclon(rec),lats(fynindex(rec)),lats(fysindex(rec)),dy(rec)

        !-------------- Linear interpolation ----------------
        !  Assign weights:

        if (same_y(rec).and.same_x(rec))then
            factx = 1.0
            facty = 1.0
        elseif (same_y(rec))then
            facty = 1.0
            factx = (frclon(rec) - lons(fxwindex(rec)))/(lons(fxeindex(rec))-lons(fxwindex(rec)))
        elseif (same_x(rec))then
            facty = (frclat(rec) - lats(fysindex(rec)))/(lats(fynindex(rec))-lats(fysindex(rec)))
            factx = 1.0
        else
            factx = (frclon(rec) - lons(fxwindex(rec)))/(lons(fxeindex(rec))-lons(fxwindex(rec)))
            facty = (frclat(rec) - lats(fysindex(rec)))/(lats(fynindex(rec))-lats(fysindex(rec)))
        endif

        print *,'frclon,lon_w,lon_e = ',frclon(rec),lons(fxwindex(rec)),lons(fxeindex(rec))
        print *,'frclat,lat_w,lat_e = ',frclat(rec),lats(fysindex(rec)),lats(fynindex(rec))

        print *,'****************'
        print *,'rec = ',rec
        print *,'****************'


        ! ---------------- Linear interpolation -------------
        ! First, create the '2d' arrays - that is, the mslp
        ! value interpolated to the user defined lat and lon.


        do i=1,maxfiles
            !	  do i=3,5

            if (i.eq.1)then

                do k=1,NMSL

                    n_2d(rec) = (dat1(fxwindex(rec),fynindex(rec),k,rec)*(1-factx)) &
                        &                    +(dat1(fxeindex(rec),fynindex(rec),k,rec)*factx)

                    s_2d(rec) = (dat1(fxwindex(rec),fysindex(rec),k,rec)*(1-factx)) &
                        &                    +(dat1(fxeindex(rec),fysindex(rec),k,rec)*factx)

                    z1_2d(rec) = (s_2d(rec)*(1-facty))+(n_2d(rec)*facty)

                    e_2d(rec) = (dat1(fxeindex(rec),fysindex(rec),k,rec)*(1-facty)) &
                        &                    +(dat1(fxeindex(rec),fynindex(rec),k,rec)*facty)

                    w_2d(rec) = (dat1(fxwindex(rec),fysindex(rec),k,rec)*(1-facty)) &
                        &                    +(dat1(fxwindex(rec),fynindex(rec),k,rec)*facty)

                    z2_2d(rec) = (w_2d(rec)*(1-factx))+(e_2d(rec)*factx)

                    i_int1(rec) = (z1_2d(rec) + z2_2d(rec))/2

                enddo

            else			! not mslp data, so assign other variables
                ! to 'var_in' working array...

                if (i.eq.2)then
                    var_in = dat2
                elseif (i.eq.3)then
                    var_in = dat3
                elseif (i.eq.4)then
                    var_in = dat4
                elseif (i.eq.5)then
                    var_in = dat5
                elseif (i.eq.6)then
                    var_in = dat6
                elseif (i.eq.7)then                              ! vjb 2.3
                    var_in = dat7
                endif


                do k=1,nlev		! ...and do the interpolation...

                    n_var(k,rec) = (var_in(fxwindex(rec),fynindex(rec),k,rec)*(1-factx)) &
                        &                    +(var_in(fxeindex(rec),fynindex(rec),k,rec)*factx)

                    s_var(k,rec) = (var_in(fxwindex(rec),fysindex(rec),k,rec)*(1-factx)) &
                        &                    +(var_in(fxeindex(rec),fysindex(rec),k,rec)*factx)

                    z1_var(k,rec) = (s_var(k,rec)*(1-facty))+(n_var(k,rec)*facty)

                    e_var(k,rec) = (var_in(fxeindex(rec),fysindex(rec),k,rec)*(1-facty)) &
                        &                    +(var_in(fxeindex(rec),fynindex(rec),k,rec)*facty)

                    w_var(k,rec) = (var_in(fxwindex(rec),fysindex(rec),k,rec)*(1-facty)) &
                        &                    +(var_in(fxwindex(rec),fynindex(rec),k,rec)*facty)

                    z2_var(k,rec) = (w_var(k,rec)*(1-factx))+(e_var(k,rec)*factx)

                    i_var(k,rec) = (z1_var(k,rec) + z2_var(k,rec))/2


                    if (i.eq.2)then
                        i_int2(k,rec) = i_var(k,rec)		    ! 2 = Z

                    elseif (i.eq.3)then
                        i_int3(k,rec) = i_var(k,rec)		    ! 3 = U

                    elseif (i.eq.4)then
                        i_int4(k,rec) = i_var(k,rec)		    ! 4 = V

                    elseif (i.eq.5)then

                        if (dx(rec).eq.0.)then
                            gradtx(k,rec) = 0.0
                        else
                            gradtx(k,rec) = (e_var(k,rec)-w_var(k,rec))/dx(rec)    ! 5 = Temp.
                        endif
                        if (dy(rec).eq.0.)then
                            gradty(k,rec) = 0.0
                        else
                            gradty(k,rec) = (n_var(k,rec)-s_var(k,rec))/dy(rec)
                        endif

                        i_int5(k,rec) = i_var(k,rec)

                    elseif (i.eq.6)then

                        if (dx(rec).eq.0.)then
                            gradqx(k,rec) = 0.0
                        else
                            gradqx(k,rec) = (e_var(k,rec)-w_var(k,rec))/dx(rec)    ! 6 = Spec.H
                        endif
                        if (dy(rec).eq.0.)then
                            gradqy(k,rec) = 0.0
                        else
                            gradqy(k,rec) = (n_var(k,rec)-s_var(k,rec))/dy(rec)
                        endif

                        i_int6(k,rec) = i_var(k,rec)

                    elseif (i.eq.7)then
                        i_int7(k,rec) = i_var(k,rec)                  ! 7 = omega    ! vjb 2.3

                    endif

                enddo

            endif

          !	  print *,'finished var(',i,')'

        enddo

    enddo

    !print *,'u,v,w = ',i_int3,i_int4,i_int7
    do i=1,nlev
        print *,'u=i_int3 for ilev=',i,nlev,nrec
        write(6,'(10(E15.6,","))')(i_int3(i,j),j=1,nrec)
    enddo


    !----- Insert routine to derive vertical advective tendancy
    !
    call gradz(nrec,nlev,levs,i_int5,i_int6,i_int7,gradtz,gradqz)

    !----- Combine X and Y (and Z) components for advective tendencies -----
    !-----------------------------------------------------------------------
    !      If vertical windspeed is NaN, then just use horizontal
    !      component
    !-----------------------------------------------------------------------
    !
    do j=1,nrec
        do i=1,nlev

            !           gradt(i,j) = -1*((i_int3(i,j)*gradtx(i,j)) + (i_int4(i,j)*gradty(i,j)))

            !	   gradq(i,j) = -1*((i_int3(i,j)*gradqx(i,j)) + (i_int4(i,j)*gradqy(i,j)))


            !-------- not used atm ---
            !
            gradt(i,j) = -1*((i_int3(i,j)*gradtx(i,j)) + (i_int4(i,j)*gradty(i,j)) + &
                &                     (i_int7(i,j)*gradtz(i,j)))
            gradq(i,j) = -1*((i_int3(i,j)*gradqx(i,j)) + (i_int4(i,j)*gradqy(i,j)) + &
                &                     (i_int7(i,j)*gradqz(i,j)))
          !
          !	if ((i_int3(i,j)*gradtx(i,j)).lt.(i_int7(i,j)*gradtz(i,j)))then
          !	 print *,'u,w,gradtz =',i_int3(i,j),gradtx(i,j),i_int7(i,j),gradtz(i,j)
          !	endif
          !	if ((i_int3(i,j)*gradqx(i,j)).lt.(i_int7(i,j)*gradqz(i,j)))then
          !	 print *,'u,w,gradqz =',i_int3(i,j),gradqx(i,j),i_int7(i,j),gradqz(i,j)
          !	endif
          !-------------

        enddo
    enddo


    !----- Tidy up and checking -----
    !
    deallocate(n_2d)			! north
    deallocate(s_2d)			! south
    deallocate(e_2d)			! east
    deallocate(w_2d)			! west
    deallocate(z1_2d)			! interpolated N/S points
    deallocate(z2_2d)			! interpolated E/W points

    deallocate(z1_var)		! interpolated N/S points
    deallocate(z2_var)		! interpolated E/W points
    deallocate(n_var)		! north
    deallocate(s_var)		! south
    deallocate(e_var)		! east
    deallocate(w_var)		! west
    deallocate(i_var)		! final interpolated output

    return

end subroutine interp_linear
