program genesis

  !  This is the genesis_main program
  !
  !  genesis intends to be a more generic rewrite of the nc2scm, scm2scum,
  !  and cre8nlist programme flow to create UM-SCM namelists specified
  !  at a user-defined location (lat/lon) taking data from NetCDF
  !  input files.
  !
  !  Parameters read from model/reanalyses are (in order): MSLP,Z,U,V,T,Q
  !  although these could be changed in later versions.  Genesis will then
  !  determine weighting factors and the data is interpolated first in the
  !  north-south direction, then the east-west direction.  From this,
  !  advective tendencies are also calculated.
  !
  !  The data for the user-defined location is then interpolated on to
  !  UM theta/rho levels and a namelist is generated from a standard
  !  template.  It is recommended the user check the namelist settings
  !  before running the SCM.
  !  - vjb 29/1/09
  !
  !       Upgrade to use generic NetCDF read routines (from diurnal2d)
  !                                       - vjb 9/11/2010


  use global, only: max_nrecs, maxfiles, num, secday, umlev, zero
  use netcdf_type

  implicit none

  !  implicit none

  !------------------------------ Declarations -------------------------
  !  Command line options
  !
  integer :: ngtopt
  integer :: iopt,nopt
  integer  :: narg
  integer  :: iargc
  character*(50)  :: optarg
  logical  :: debug,date,offset
  logical  :: rh
  logical  :: pascal
  character*25  :: infile(num)

  logical  :: there
  logical  :: L_vlatlon
  logical  :: L_w_inOK


  !------------------------------
  !  More command line parameters
  !
  !  character*4  :: chnt  ! time increments
  !  character*2  :: chnlev  ! # levels of input data
  !  character*5  :: chx    ! # lons, user longitude !glr
  character*7  :: usrx    ! # lons, user longitude !glr
  !  character*5  :: chy    ! # lats, user latitude !glr
  character*7  :: usry          ! # lats, user latitude !glr
  integer  :: nfiles  ! # input nc files
  integer  :: ndates  ! # dates in datefile
  integer,dimension(:),allocatable  :: ymd,hr  ! date components

  logical,dimension(:),allocatable  :: same_x,same_y   ! interp check
  !  logical  :: offset
  !  logical  :: date

  real  :: usrlat,usrlon ! user lat/lon

  !------------------------------
  !  Counters
  !
  integer  :: i,j,k,crec

  !-----------------------------
  !  NetCDF settings
  !

  type(netcdf_metadata) :: nc_meta
  integer  :: NRECS,NLVLS,NLONS,NLATS

  character*8  :: chymd(max_nrecs)  ! user def. date

  !-----------------------------
  !  Data arrays

  real,dimension(:,:,:,:),allocatable  :: data
  

  real,dimension(:),allocatable  :: rec

  real,dimension(:,:,:,:),allocatable  :: msl_in
  real,dimension(:,:,:,:),allocatable  :: z_in
  real,dimension(:,:,:,:),allocatable  :: u_in
  real,dimension(:,:,:,:),allocatable  :: v_in
  real,dimension(:,:,:,:),allocatable  :: t_in
  real,dimension(:,:,:,:),allocatable  :: q_in
  real,dimension(:,:,:,:),allocatable  :: w_in                   ! vjb 2.3

  real,dimension(:),allocatable  :: msl
  real,dimension(:,:),allocatable  :: z
  real,dimension(:,:),allocatable  :: u
  real,dimension(:,:),allocatable  :: v
  real,dimension(:,:),allocatable  :: t
  real,dimension(:,:),allocatable  :: q
  real,dimension(:,:),allocatable  :: w                      ! vjb 2.3

  real,dimension(:),allocatable  :: lats
  real,dimension(:),allocatable  :: lons
  real,dimension(:),allocatable  :: levs

  real  :: adum

  !-----------------------------
  !  Variable lon/lat forcing arrays

  real,dimension(:,:),allocatable  :: gradt
  real,dimension(:,:),allocatable  :: gradq

  real,dimension(:,:),allocatable  :: ug
  real,dimension(:,:),allocatable  :: vg
  real,dimension(:,:),allocatable  :: dug
  real,dimension(:,:),allocatable  :: dvg
  real  :: f
  real,dimension(:,:),allocatable  :: pt

  real,dimension(:,:),allocatable  :: p_um,pt_um
  real,dimension(:,:),allocatable  :: u_um,v_um
  real,dimension(:,:),allocatable  :: ug_um,vg_um
  real,dimension(:,:),allocatable  :: t_um,q_um
  real,dimension(:,:),allocatable  :: w_um                  ! vjb 2.3
  real,dimension(:,:),allocatable  :: gradt_um,gradq_um

  real,dimension(:),allocatable  :: frclat
  real,dimension(:),allocatable  :: frclon
  real,dimension(:),allocatable  :: dx,dy
  integer,dimension(:),allocatable  :: frcynindex
  integer,dimension(:),allocatable  :: frcysindex
  integer,dimension(:),allocatable  :: frcxeindex
  integer,dimension(:),allocatable  :: frcxwindex
  integer  :: nfrc
  integer, dimension(num) :: varid
  integer, dimension(4, num) :: var_dimids


  !-------------------- Command line setup ------------------
  !  Usage options
  narg = iargc()

  if (narg.eq.0)then
    call write_help()
    stop
  endif

  !----------------------------------------------
  !  Set defaults and command line options

  debug = .false.
  date = .false.
  rh = .false.
  pascal = .false.
  !  chnt = '150'
  !  chnlev = '20'
  !  chx = '29'
  !  chy = '17'
  offset = .false.
  !  same_y = .false.
  !  same_x = .false.
  usrx = 'A'
  usry = 'B'

  do while (iopt.ne.-1)
    iopt = ngtopt('dhX:Y:RMDO',nopt,optarg)  ! <<
    if (char(iopt).eq.'h') then
      call write_help()
      stop
    end if
    if (char(iopt).eq.'d') debug = .true.
    if (char(iopt).eq.'X') usrx = trim(optarg)
    if (char(iopt).eq.'Y') usry = trim(optarg)
    if (char(iopt).eq.'R') rh = .true.
    if (char(iopt).eq.'M') pascal = .true.
    if (char(iopt).eq.'D') date = .true.
    if (char(iopt).eq.'O') offset = .true.
  enddo

  if (debug) print *,'Debug option selected'
  if (date) print *,'User defined date file specified'
  if (rh) print *,'Convert RelH to specific humidity flag raised'
  if (pascal) print *,'Convert surface pressure from hPa to Pa flag raised'
  if (offset) print *,'User defined offset file specified'

  !------ Interrogate filesfile -------
  !
  inquire(file='files.inp',exist=there)
  if(.not.there) stop' files.inp not found'

  inquire(file='base.inp',exist=there)
  if(.not.there) stop' base.inp namelist file not found'

  inquire(file='template.scm',exist=there)
  if(.not.there) stop'ERROR: template.scm file not found in run directory'

  open (1,file='files.inp')

  do i=1,num
    read (1,*,end=9)infile(i)
    inquire(file=infile(i),exist=there)
    if(.not.there) stop'input NetCDF not found'
    if (debug) print *,'i,infile(i)=',i,infile(i)
  enddo

9 continue

  nfiles = i-1

  L_w_inOK   = .false. !if we do/do not have w_in,set w_inOK=true/false
  if (debug) print *,'L_w_inOK = ',L_w_inOK
  if (debug) print *,'nfiles,maxfiles = ',nfiles,maxfiles
  if (nfiles.lt.maxfiles.and.(L_w_inOK)) stop'Not enough input files specified'
  if (nfiles.gt.maxfiles) stop'Too many input files specified'


  !----- Set up command line parameters (if req'd)
  !

  if (usrx.eq.'A') stop'User longitude not specified (eg. 243.7)'
  if (usry.eq.'B') stop'User latitude not specified (eg. 65.0)'

  read (usry,'(f7.3)')usrlat                                          !glr change from 5.1

  if (usrlat.lt.-90.or.usrlat.gt.90)stop'ERROR: Latitude outside &
    &       geographical boundaries'

  if (debug) print *,'usrlat = ',usrlat

  read (usrx,'(f7.3)')usrlon                                          !glr change from 5.1

  if (usrlon.lt.0.or.usrlon.gt.360)stop'ERROR: Longitude outside &
    &       geographical boundaries'

  if (debug) print *,'usrlon = ',usrlon

  !------- Open NetCDFs -------------------
  !

  if (debug) print *,'Calling ncread_dim...'

  do k=1,nfiles

    if (k.gt.1)then
      deallocate(lons)
      deallocate(lats)
      deallocate(levs)
      deallocate(rec)
      deallocate(data)
    endif

    call ncread_dim(nc_meta,debug,infile(k), varid(k), &
    & var_dimids(:, k))

    NLONS = nc_meta%lon_n
    NLATS = nc_meta%lat_n
    NLVLS = nc_meta%lvl_n
    NRECS = nc_meta%rec_n

    allocate(lons(NLONS))
    allocate(lats(NLATS))
    allocate(levs(NLVLS))
    allocate(rec(NRECS))
    allocate(data(NLONS,NLATS,NLVLS,NRECS))

    if (debug) print *,'nlon,nlat,nlev,nrec = ',NLONS,NLATS,NLVLS,NRECS

    call ncread_data(nc_meta,lons,lats,levs,rec,data, &
    &   var_dimids(:, k))

    if (k.eq.1)then

      allocate(msl_in(NLONS,NLATS,NLVLS,NRECS))

      msl_in = data

      if (pascal)then
        msl_in = msl_in*100                          ! convert MSLP from hPa to Pa
      endif

      if (debug) print *,'msl = ',msl_in(1,20,1,1)

    elseif (k.eq.2)then

      allocate(z_in(NLONS,NLATS,NLVLS,NRECS))

      z_in = data

      if (debug) print *,'z = ',z_in(1,10,1,1)

    elseif (k.eq.3)then

      allocate(u_in(NLONS,NLATS,NLVLS,NRECS))

      u_in = data

      if (debug) print *,'u = ',u_in(1,1,1,1)

    elseif (k.eq.4)then

      allocate(v_in(NLONS,NLATS,NLVLS,NRECS))

      v_in = data

      if (debug) print *,'v = ',v_in(1,1,1,1)

    elseif (k.eq.5)then

      allocate(t_in(NLONS,NLATS,NLVLS,NRECS))

      t_in = data

      if (debug) print *,'t = ',t_in(1,1,1,1)

    elseif (k.eq.6)then

      allocate(q_in(NLONS,NLATS,NLVLS,NRECS))

      q_in = data

      if (debug) print *,'q = ',q_in(20,20,10,1)
      if (.not. L_w_inOK) then
        print *,'L_w_inOK',L_w_inOK
        allocate(w_in(NLONS,NLATS,NLVLS,NRECS))           ! vjb 2.3
        w_in = 0.0*q_in
        if (debug) print *,'w = ',w_in(20,20,10,1)
      endif

    elseif (k.eq.7)then

      allocate(w_in(NLONS,NLATS,NLVLS,NRECS))           ! vjb 2.3

      w_in = data

      if (debug) print *,'w = ',w_in(20,20,10,1)

    endif

  enddo


  !------- Note for later --------
  !
  ! If planning on reading in ERA fields then note that the
  ! fields are oriented differently to UM fields:
  !  lat(ERA): 90 --> -90     lat(UM): -90 --> 90
  !  lon(ERA): 0 --> 360      lon(UM): 0 --> 360
  !  lev(ERA): TOA --> 1000mb lev(UM): 1000mb --> TOA
  !
  ! Also, the MSL data for ERA only hs 3 dimensions making
  ! it awkward to read using ncread_data atm (grr..)
  ! (this has now been updated with the -R option in the command line, see help). vjb
  !-------------------------------


  print *,'returned from ncread'


  !------ If required, convert relative humidity to specific humidity here ------
  !

  if (rh)then

    print *,'-R flag raised: converting relative humidities to specific humidities..'

    call relh2q(NLONS,NLATS,NLVLS,NRECS,levs,t_in,q_in)

  endif


  !
  !------- glr added bit --------
  !------- Read the variable lon/lat locations from lonlat.dat file ------------
  !
  L_vlatlon = .false.  !assume single lonlat as in input
  !       L_vlatlon = .true.   !read in lonlat data

  if (L_vlatlon)then

    inquire(file='lonlat.dat',exist=there)
    if(.not.there) stop' lonlat.dat file not found'

    open(1,file='lonlat.dat')

    do i=1,max_nrecs  ! Next, check to make sure
      read(1,*,end=998)chymd(i)  ! the user-defined lonlat file
    enddo  ! matches the number of records

998 continue  ! in the NetCDF file

    close(1)

    nfrc = i-1

    if (debug) print *,'nfrc = ',nfrc

    !-------------------------------------------------------------------------------
    if (nfrc.ne.NRECS)stop'ERROR: user defined lonlat mismatch with NRECS'
    !-------------------------------------------------------------------------------
    !
    allocate(frclon(NRECS))  !
    allocate(frclat(NRECS))  ! allocate the lonlat arrays
    open(1,file='lonlat.dat')              ! re-initialise read of date file

    do i=1,NRECS
      read(1,'(1x,F7.3,1x,F7.3)')frclon(i),frclat(i)  !! << watch the formatting !!
    enddo

  else                                   ! variable lat/lon is false

    nfrc = 1

    allocate(frclon(NRECS))  !
    allocate(frclat(NRECS))  ! allocate the lonlat arrays

    do i=1,NRECS

      read (usrx,'(f7.3)')frclon(i)                   !! << watch the formatting !!
      read (usry,'(f7.3)')frclat(i)                    !! << watch the formatting !!

    enddo

  endif

  close(1)

  print *,'usrlon,usrlat=',usrlon,usrlat
  print *,'frclon=',frclon
  print *,'frclat=',frclat


  !------- Find nearest lat/lon interval ----------------
  !
  print *,'Finding nearest neighbour grid points'
  print *,'nfrc = ',nfrc

  allocate(frcynindex(NRECS))  !
  allocate(frcysindex(NRECS))  !
  allocate(frcxeindex(NRECS))  !
  allocate(frcxwindex(NRECS))  !

  call find_latlon(NRECS,NLATS,NLONS,frclat,frclon,lats,lons,frcynindex,frcysindex, &
    &                     frcxeindex,frcxwindex,debug)

  if (debug) write(*,*)'usrlat,usrlon = ',usrlat,usrlon


  !------- Set levs array to be in Pascals ---------------
  !
  do k=1,NLVLS
    levs(k) = 100*levs(k)
  enddo

  !------- Interpolate data to user-defined lat/lon ------
  !-------------------------------------------------------
  !----- Insert alternative interpolation method here ----
  !--------------------- NOT USED ------------------------
  !
  !-------------- Linear interpolation ----------------
  !  Now we have the nearest neighbour grid points to the
  !  user defined lat/lon ready for linear interpolation.
  !  (ynindex,ysindex,xeindex,xwindex)
  !

  allocate(msl(NRECS))
  allocate(z(NLVLS,NRECS))
  allocate(u(NLVLS,NRECS))
  allocate(v(NLVLS,NRECS))
  allocate(t(NLVLS,NRECS))
  allocate(q(NLVLS,NRECS))
  allocate(w(NLVLS,NRECS))                               ! vjb 2.3
  allocate(gradt(NLVLS,NRECS))
  allocate(gradq(NLVLS,NRECS))
  allocate(same_y(NRECS))
  allocate(same_x(NRECS))
  allocate(dx(NRECS))
  allocate(dy(NRECS))

  !        print *,'glr11a0 shape2',shape(q),shape(t),shape(gradt),shape(gradq)

  do i=1,NRECS
    same_y(i) = .false.
    if (frcynindex(i).eq.frcysindex(i))then
      same_y(i) = .true.
    endif

    same_x(i) = .false.
    if (frcxeindex(i).eq.frcxwindex(i))then
      same_x(i) = .true.
    endif
  enddo

  call interp_linear(NLATS,NLONS,NLVLS,NRECS,frcynindex,frcysindex, &
    &                     frcxeindex,frcxwindex,frclat,frclon,lats,lons,levs, &
    &                     same_x,same_y,msl_in,z_in,u_in,v_in,t_in,&
    &                     q_in,w_in,msl,z,u,v,t,q,w,dx,dy,gradt,gradq,debug)

  print *,'glr11a2 shape2',shape(p_um),shape(pt_um),shape(gradt),shape(gradq)

  !  print *,'msl,z,u,v,t,q,w =',u             ! vjb 2.3


  !------- Calculate advective tendencies ---------------
  !

  print *,'Calculating advective tendencies'
  print *,'Temperature...'

  !        print *,'glr11a3 shape2',shape(p_um),shape(pt_um),shape(gradt),shape(gradq)
  do j=1,NRECS
    do i=1,NLVLS
      gradt(i,j) = gradt(i,j)*secday  ! change units to K/day
    enddo
  enddo

  !  if (debug) print *,'gradt = ',(gradt(i,1),i=1,NLVLS)

  !  -----------------------------

  print *,'Specific humidity...'

  do j=1,NRECS
    do i=1,NLVLS
      gradq(i,j) = gradq(i,j)*secday     ! change units to kgkg/day
    enddo
  enddo

  !  if (debug) print *,'gradq = ',(gradq(i,1),i=1,NLVLS)

  !  stop

  !------- Derive geostrophic wind profiles --------------
  !
  allocate(ug(NLVLS,NRECS))
  allocate(vg(NLVLS,NRECS))
  allocate(pt(NLVLS,NRECS))

  call geostr(frclat,NRECS,NLVLS,NLATS,NLONS,levs,   &
    &        frcynindex,frcysindex,frcxeindex,frcxwindex,t,z_in,dx,dy,ug,vg,f,pt,debug)

  !  if (debug) print *,'ug = ',((ug(k,j),k=1,NLVLS),j=1,NRECS)
  !  if (debug) print *,'f = ',f


  !------- Assign dates/times to time records ------------
  !
  if (date)then

    inquire(file='dates.dat',exist=there)
    if(.not.there) stop' dates.dat file not found'

    open(1,file='dates.dat')

    do i=1,max_nrecs  ! Next, check to make sure
      read(1,*,end=999)chymd(i)  ! the user-defined dates file
    enddo  ! matches the number of records
  ! in the NetCDF file - and
999 continue  ! fill chymd whilst we're here...

    close(1)

    ndates = i-1

    if (debug) print *,'ndates = ',ndates
    !-------------------------------------------------------------------------------
    if (ndates.ne.NRECS)stop'ERROR: user defined dates mismatch with NRECS'
    !-------------------------------------------------------------------------------

    allocate(ymd(NRECS))  !
    allocate(hr(NRECS))  ! allocate the date arrays

    open(1,file='dates.dat')  ! re-initialise read of date file

    do i=1,NRECS
      read(1,'(1x,I8,2x,I4)')ymd(i),hr(i)  !! << watch the formatting !!
    enddo

  else

    print *,'**********'
    stop'ERROR: Dates must be assigned to input fields using dates.dat file (-D option)'

    allocate(ymd(NRECS))  !
    allocate(hr(NRECS))  ! allocate the date arrays

    do i=1,NRECS  ! NOTE: For the time increment calculation this
      ymd(i) = 0  ! must be set to numbers that are not zero
      hr(i) = 0  ! therefore this option has been deactivated
    enddo

  endif


  close(1)

  !  write(*,*)"z = ",((z(k,crec),k=1,NLVLS),crec=1,NRECS)



  !------------ Write out *.scm files ----------
  !
  print *,'Writing *.scm file'

  open(3,file='genesis.scm')

  do crec=1,NRECS

    if (date)then
      write(3,*)ymd(crec),hr(crec),msl(crec),NLVLS+1,usrlat,usrlon
    else
      write(3,*)crec,crec,msl(crec),NLVLS+1,usrlat,usrlon           ! deactivated
    endif

    write(3,20)zero,msl(crec),t(1,crec),pt(1,crec), &  ! sfc data
      &        q(1,crec),u(1,crec),v(1,crec),w(1,crec),ug(1,crec),vg(1,crec)

    do k=1,NLVLS                      !profile data on input netcdf levels NLVLS
      write(3,20)z(k,crec),levs(k),t(k,crec),pt(k,crec), &
        &         q(k,crec),u(k,crec),v(k,crec),w(k,crec),ug(k,crec),vg(k,crec)
    enddo

  enddo

20 format(2(f9.0,' '),2(f8.2),' ',e12.5,' ',5(f8.2))

  close(3)

  print *,'Outfile written successfully: genesis.scm'

  print *,'Searching for base namelist...did you specify it in cmd line?'
  crec=NRECS !glr 29-04-2014


  !--------------------------------------------------------------------
  !------ Interpolate profiles to Charney-Phillips vertical grid ------
  !== interpolate from input netcdf NLVLS to output um umlevs

  allocate(p_um(umlev+1,NRECS))
  allocate(pt_um(umlev+1,NRECS))
  allocate(u_um(umlev,NRECS))
  allocate(ug_um(umlev,NRECS))
  allocate(v_um(umlev,NRECS))
  allocate(vg_um(umlev,NRECS))
  allocate(w_um(umlev+1,NRECS))
  allocate(t_um(umlev+1,NRECS))
  allocate(q_um(umlev+1,NRECS))
  allocate(gradt_um(umlev+1,NRECS))
  allocate(gradq_um(umlev+1,NRECS))

  !       setup initial um field values
  do j=1,NRECS
    do k=1,umlev
      ug_um(k,j)=0.0
      vg_um(k,j)=0.0
      u_um(k,j) =0.0
      v_um(k,j) =0.0
    enddo
    do k=1,umlev+1
      p_um(k,j) =0.0
      pt_um(k,j) =0.0
      w_um(k,j) =0.0
      t_um(k,j) =0.0
      q_um(k,j) =0.0
      gradt_um(k,j) =0.0
      gradq_um(k,j) =0.0
    enddo
  enddo

  call charney_phillips(NRECS,NLVLS,levs,z,u,v,w,t,q,pt,msl,gradt,gradq, &
    &      u_um,v_um,w_um,t_um,q_um,p_um,pt_um, &
    &                              gradt_um,gradq_um)

  !  call charney_phillips(NRECS,NLVLS,levs,z,u,ug,v,vg,w,t,q,pt,msl,gradt,gradq, &
    !     &      u_um,ug_um,v_um,vg_um,w_um,t_um,q_um,p_um,pt_um, &
    !     &                              gradt_um,gradq_um,debug)

  open(4,file='charney.scm')

  do k=1,umlev+1
    adum=0.0
    if (k.le.NLVLS) adum=levs(k)
    write(4,21)p_um(k,crec),adum,t_um(k,crec),pt_um(k,crec), &
      &         q_um(k,crec),u_um(k,crec),v_um(k,crec),w_um(k,crec)
    print *,'D03',k,crec,v_um(k,crec),w_um(k,crec)
  enddo

  !21      format(1(f9.0,' '),2(f8.2),' ',e12.5,' ',3(f8.2))
21 format(1(f9.0,' '),3(f9.2),' ',e12.5,' ',3(f8.2))
  close(4)

  print *,'Outfile written successfully: charney.scm'

  !------ Create UM-SCM namelist ------
  !

  print *,'Creating output namelist...'

  !        print *,'glr11a shape1',shape(u_um),shape(v_um),shape(t_um),shape(q_um)
  !        print *,'glr11a shape2',shape(p_um),shape(pt_um)
  !        print *,'glr11a shape3',shape(gradt),shape(gradq),shape(gradt_um),shape(gradq_um)
  !        print *,'glr11a gen NRECS,umlev ',NRECS,umlev
  !        print *,'glr11a gen u_um ',u_um
  !        print *,'glr11a gen v_um ',v_um
  !        print *,'glr11a gen t_um ',t_um
  !        print *,'glr11a gen q_um ',q_um
  !        print *,'glr11a gen p_um ',p_um
  !        print *,'glr11a gen pt_um ',pt_um
  !        print *,'glr11a gen gradt_um ',gradt_um
  !        print *,'glr11a gen gradq_um ',gradq_um
  !        print *,'glr11a gen ymd,hr,usrlat,usrlon,debug ',ymd,hr,usrlat,usrlon,debug

  !------- Calculate U and V geostrophic advective tendencies ------
  !
  allocate(dug(NLVLS,NRECS))
  allocate(dvg(NLVLS,NRECS))

  do j=1,NRECS
    do k=1,NLVLS

      dug(k,j) = f * (vg_um(k,j) - v_um(k,j))
      dvg(k,j) = f * (ug_um(k,j) - u_um(k,j))

    enddo
  enddo

  if (debug) print *,'dug,dvg = ',dug(10,2),dvg(10,2)


  call create_namelist(NRECS,umlev,u_um,dug,v_um,dvg,w_um,t_um,q_um,p_um,pt_um,  &
    &       gradt_um,gradq_um,ymd,hr,usrlat,usrlon,debug)

contains
  subroutine write_help()
    implicit none
    write(*,*)''
    write(*,*)'Usage: genesis'
    write(*,*)'  genesis [options] < base.inp'
    write(*,*)'  -d        verbose, debug mode'
    write(*,*)'  -h        help'
    write(*,*)'  -X        user defined longitude (eg: 243.7)'
    write(*,*)'  -Y        user defined latitude (eg: 65.0)'
    write(*,*)'  -R        convert RelH to spec hum (def = false)'
    write(*,*)'  -M        convert sfc pressure from hPa to Pa (def = false)'
    write(*,*)'  -D        use date file (def = false)'
    write(*,*)'  -O        use offset file (def = false)'
    write(*,*)'  -x        number of longitude values in netcdf file'
    write(*,*)'  -y        number of latitude values in netcdf file'
    write(*,*)'  -l        number of vertical levels in netcdf file'
    write(*,*)'  -t        number of timesteps in netcdf file'
    write(*,*)'  -U < file add namelists for various parameters.'
    write(*,*)''
    write(*,*)'If the only moisture field available is Relative Humidity, the -R option'
    write(*,*)'will calculate the equivalent specific humidity from the temperature field'
    write(*,*)''
    write(*,*)'You also need a few files in the run directory for this to work.'
    write(*,*)' a) files.inp - a file listing the *.nc file names in order: ps,z,u,v,t,q'
    write(*,*)'                one filename per line'
    write(*,*)' b) offset.dat (optional) - list of offset parameters from'
    write(*,*)'     ncdump hdr [scale, add]'
    write(*,*)'              - if geopotential needs converting only (m2/s2 >> m)'
    write(*,*)'                then use -O with offset file set to 1 0 for each.'
    write(*,*)' c) dates.dat (reqd) - a list of dates/times' // &
      &                             ' corresponding to data records.'
    write(*,*)'     [generate using datetraj - yyyymmdd hhhh]'
    write(*,*)''
    write(*,*)'NB: ensure all fields are on the same grids - check for staggering of'
    write(*,*)'    momentum/scalar fields.  Recommend Xconv to do quick conversion.'
    write(*,*)''
    write(*,*)' Example: '
    write(*,*)'  ./genesis -d -t109 -l37 -x27 -y27 -X148.2 -Y-35.7 -D -U < base.inp'
  end subroutine write_help
end program genesis
