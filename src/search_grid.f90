subroutine search_grid(usrpt,dat,max_array_size,positive,grid_index)

  !	Given a grid point (usrpt), will find nearest neighbour in dat array.
  !	Choose positive = .true. to search 'up' and .false. to search 'down'.

  implicit none

  !------- declarations

  real,parameter  :: acc = 0.007 ! accuracy of real number search
  integer         :: max_array_size ! size of destination array
  real,intent(in) :: usrpt ! user defined grid point (f5.1)
  real,dimension(max_array_size),intent(in)  :: dat ! destination grid array
  logical     :: positive ! search up or down (+/-)
  integer,intent(out)             :: grid_index ! output grid index
  real        :: dd,delta,degsearch,gacc
  integer     :: mmax

  !------- working

  integer     :: m,n    ! counters
  real        :: gridpt ! working grid point value

  !---------------- Search -------------
  !

  degsearch = 5.0 !valid for access-g
  degsearch = 0.5 !valid for access-c where dlat=0.05deg =>steps of 0.001
  degsearch = 5.0 !valid for era      where dlat=0.75deg =>steps of 0.01
  mmax = 500
  delta = degsearch/float(mmax)
  gacc = 0.5*delta
  if (positive)then

    gridpt = usrpt        ! initialise
    up_loop: do m=1,mmax           ! search 'up' 5 deg in steps of 0.01
      do n=1,max_array_size
        dd = gridpt-dat(n)
        if (abs(dd).le.gacc)then
          exit up_loop
        endif
      enddo
      gridpt = gridpt+delta
    enddo up_loop

    grid_index = n

  else      ! search 'down' 5 deg in steps of 0.01

    gridpt = usrpt
    down_loop: do m=1,mmax
      do n=max_array_size,1,-1
        dd = gridpt-dat(n)
        if (abs(dd).le.gacc)then
          exit down_loop
        endif
      enddo
      gridpt = gridpt-delta
    enddo down_loop

    grid_index = n

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
