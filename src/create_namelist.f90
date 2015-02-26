
subroutine create_namelist(nrec,nlev,u_um,dug,v_um,dvg,w_um,t_um,q_um, &
    &				   p_um,pt_um,gradt_um,gradq_um,ymd,hr, &
    &                             usrlat,usrlon,debug)

  !	Reads in user-defined settings from input namelist and creates
  !	an output namelist based upon generic GABLS settings
  !								- vjb 8/5/09
  !
  !	nb: - tstar uses a default of 288K, must be user defined by direct edit
  !           - wi and w_inc also not extracted from NWP, also user defined
  !           - code may need modifying if surface forcing by E/H fluxes desired
  !	    - user can add own fields, copy lines labelled 21-23 and edit
  !	      and include an extra namelist logic switch if desired
  !
  use global

  implicit none

  !------ Input variables ------
  !
  integer,intent(in)				:: nrec,nlev
  real,dimension(umlev,nrec),intent(in)		:: u_um,v_um
  real,dimension(umlev,nrec),intent(in)		:: dug,dvg
  real,dimension(umlev+1,nrec),intent(in)		:: w_um,t_um,q_um,pt_um
  real,dimension(umlev,nrec),intent(in)		:: p_um
  real,dimension(umlev+1,nrec),intent(in)		:: gradt_um,gradq_um
  integer,dimension(nrec),intent(in)		:: ymd,hr
  real,intent(in)					:: usrlat,usrlon
  logical,intent(in)				:: debug

  !------ Working parameters and arrays ------
  !
  integer				:: i,j,k,rhs,t0
  integer				:: jt0,jt1,it,jt
  integer				:: nfor,nn
  real				:: dt,fact
  integer,parameter		:: bmax=10000
  character*85			:: base
  character*85			:: bch90
  character*7			:: ch7
  character*1			:: a
  integer				:: nstr,ymd1
  character*20			:: bfld
  logical				:: field,diag

  !------ Namelist variables ------
  !------ Base ------
  !
  integer                         :: sdate,shour
  integer                         :: edate,ehour
  integer				:: year,month,day

  logical				:: ui,vi,wi
  logical				:: theta,qi,p_in

  logical				:: L_windrlx
  logical				:: tau_rlx
  logical				:: L_vertadv
  logical				:: tstar_forcing
  logical				:: flux_e,flux_h
  logical				:: u_inc,v_inc,w_inc
  logical				:: t_inc,q_star
  logical				:: ichgf

  character*100			:: namelist_template

  !------ Read namelists -------
  !

  namelist/usrfields_1/ui,vi,wi,theta,qi,p_in
  namelist/usrfields_2/L_windrlx,tau_rlx,L_vertadv,tstar_forcing, &
    &                       flux_e,flux_h,u_inc,v_inc,w_inc,t_inc,     &
    &			     q_star,ichgf
  !        namelist/usrfields_3/namelist_template
  namelist/time/sdate,shour,edate,ehour,year,month,day

  read(*,usrfields_1)
  read(*,usrfields_2)
  !	read(*,usrfields_3)
  read(*,time)

  !	if (L_windrlx) print *,'L_windrlx flag raised.'
  !	if (u_inc) print *,'u_inc flag raised.'
  !	if (v_inc) print *,'v_inc flag raised.'
  !	print *,namelist_template,sdate,shour,edate,ehour


  !------  Set length of records, time intervals, etc... ------
  !
  !        print *,'glr11b shape1',shape(u_um),shape(v_um),shape(t_um),shape(q_um)
  !        print *,'glr11b shape2',shape(p_um),shape(pt_um),shape(gradt_um),shape(gradq_um)
  !        print *,'glr11b shape3',shape(gradt_um),shape(gradq_um)
  !        print *,'glr11b gen nrec,nlev ',nrec,nlev
  !        print *,'glr11b gen u_um ',u_um
  !        print *,'glr11b gen v_um ',v_um
  !        print *,'glr11b gen t_um ',t_um
  !        print *,'glr11b gen q_um ',q_um
  !        print *,'glr11b gen p_um ',p_um
  !        print *,'glr11b gen pt_um ',pt_um
  !        print *,'glr11b gen gradt_um ',gradt_um
  !        print *,'glr11b gen gradq_um ',gradq_um

  !        print *,'glr11b gen ymd,hr,usrlat,usrlon,debug ',ymd,hr,usrlat,usrlon,debug

  jt0 = 0			! period start/end indices
  jt1 = 0
  !        print *,'glr: nrec,sdate,shour,edate,ehour=',nrec,sdate,shour,edate,ehour

  do it=1,nrec
    if (jt0.eq.0.and.ymd(it).eq.sdate.and.hr(it).ge.shour)then
      jt0 = it
    elseif (jt0.eq.0.and.ymd(it).gt.sdate)then
      jt0 = it
    elseif(jt1.eq.0.and.ymd(it).eq.edate.and.hr(it).ge.ehour)then
      jt1 = it
    elseif (jt1.eq.0.and.ymd(it).gt.edate)then
      jt1 = it
    endif
    !         print *,'glr: ymd(it),hr(it),jt0,jt1=',ymd(it),hr(it),jt0,jt1
  enddo
  if (debug) print *,'jt0,jt1 = ',jt0,jt1

  dt = hr(jt0+1) - hr(jt0)
  !        print *,'glr: dt=',dt
  !        dt = dt*0.01 !glr
  !        print *,'glr: dt=dt*.01=',dt
  if (dt.le.0) dt=dt+24.0
  fact = 24.0/dt
  !	if (L_windrlx)then
  nfor = 1+jt1-jt0
  !	else
  !	 nfor = jt1-jt0
  !	endif
  nn = nfor*nlev
  if (debug) print *,'dt,fact = ',dt,fact


  !------  Read basis file, find variable name by searching for '=' ------
  !
  open(1,file='template.scm')
  open(2,file='namelist.scm')
  diag = .false.
  field=.false.

  do i=1,bmax

    read(1,'(A)',end=9)base		! read whole line
    !!
    !         print *,'base=',base
    !         print *,'shape(base)=',shape(base)
    !         print *,'len(base)=',len(base)
    !         print *,'base=',base
    !         print *,'len(base)=',len(base)
    !         print *,'trim(base)=',trim(base)
    !         print *,'len(trim(base))=',len(trim(base))
    !!
    backspace(1)
    read(1,*,end=9)bch90			! free format, read first string
    do j=len(bch90),1,-1
      a = bch90(j:j)
      if (a.eq.'=')then
        field=.true.
        go to 50	                        ! found field >> name it 'bfld'
      endif
    enddo
    diag = .true.
    go to 89

    50       continue

    nstr = j-1
    bfld = bch90(1:nstr)

    if (L_windrlx)then				! dump profile output

      print *,'Writing outputs as profiles'

      if (bfld.eq.'L_windrlx'.and.L_windrlx)then
        write(2,*)'L_windrlx= .TRUE.'
        field = .false.
      elseif (bfld.eq.'nfor')then
        write(2,*)'nfor=     ',1+jt1-jt0,','
        field = .false.
      elseif (bfld.eq.'TAU_RLX'.and.tau_rlx)then
        write(2,'("TAU_RLX = ",f8.0)')dt*3600.0
        field = .false.
      elseif (bfld.eq.'U_INC'.and.u_inc)then
        write(2,*)'U_INC='
        write(2,'(10(E15.6,","))')((u_um(k,jt),jt=jt0,jt1),k=1,nlev)
        field = .false.
      elseif (bfld.eq.'V_INC'.and.v_inc)then
        write(2,*)'V_INC='
        write(2,'(10(E15.6,","))')((v_um(k,jt),jt=jt0,jt1),k=1,nlev)
        field = .false.

      elseif (bfld.eq.'W_INC'.and.w_inc)then 	! vjb 2.4
        write(2,*)'W_INC='
        write(2,'(10(E15.6,","))')((w_um(k,jt),jt=jt0,jt1),k=1,nlev)
        field = .false.

      elseif (bfld.eq.'T_INC'.and.t_inc)then
        write(2,*)'T_INC='
        write(2,'(10(E15.6,","))')(((t_um(k,jt+1) - t_um(k,jt))*fact, &
          &                                 jt=jt0,jt1-1),k=1,nlev+1)
        field = .false.
      elseif (bfld.eq.'Q_STAR'.and.q_star)then
        write(2,*)'Q_STAR='
        write(2,'(10(E15.6,","))')(((q_um(k,jt+1) - q_um(k,jt))*fact, &
          &                                 jt=jt0,jt1-1),k=1,nlev+1)
        !   print *,'glr create gradq_um ',L_windrlx,((gradq_um(k,jt),jt=jt1-1,jt1-1),k=1,nlev+1)
        field = .false.

        !          elseif (bfld.eq.'T_INC'.and.t_inc)then
        !           write(2,*)'T_INC='
        !           write(2,'(10(E15.6,","))')((gradt_um(k,jt),jt=jt0,jt1),k=1,nlev+1)
        !           field = .false.
        !	  elseif (bfld.eq.'Q_STAR'.and.q_star)then
        !           write(2,*)'Q_STAR='
        !           write(2,'(10(E15.6,","))')((gradq_um(k,jt),jt=jt0,jt1),k=1,nlev+1)
        !   print *,'glr create gradq_um ',L_windrlx,((gradq_um(k,jt),jt=jt1-1,jt1-1),k=1,nlev+1)
        !           field = .false.


      endif

    else						! dump increments

      if (bfld.eq.'L_windrlx')then
        print *,'L_windrlx = .false. - writing increments/advective tendencies'
        write(2,*)'L_windrlx= .FALSE.'
        field = .false.

      elseif (bfld.eq.'nfor')then
        write(2,*)'nfor=     ',1+jt1-jt0,','         ! vjb 2.4 fix for SCM
        field = .false.

      elseif (bfld.eq.'U_INC'.and.u_inc)then
        write(2,*)'U_INC='
        write(2,'(10(E15.6,","))')(((u_um(k,jt+1) - u_um(k,jt) - dvg(k,jt))*fact, &
          &                                  jt=jt0,jt1-1),k=1,nlev)            ! - dvg(k,jt))
        field = .false.
      elseif (bfld.eq.'V_INC'.and.v_inc)then
        write(2,*)'V_INC='
        write(2,'(10(E15.6,","))')(((v_um(k,jt+1) - v_um(k,jt) + dug(k,jt))*fact, &
          &                                  jt=jt0,jt1-1),k=1,nlev)            ! + dug(k,jt))
        field = .false.

      elseif (bfld.eq.'W_INC'.and.w_inc)then 	! vjb 2.4
        write(2,*)'W_INC='
        write(2,'(10(E15.6,","))')(((w_um(k,jt+1) - w_um(k,jt))*fact, &
          &                                  jt=jt0,jt1),k=1,nlev)
        field = .false.

      elseif (bfld.eq.'T_INC'.and.t_inc)then
        write(2,*)'T_INC='
        write(2,'(10(E15.6,","))')((gradt_um(k,jt),jt=jt0,jt1),k=1,nlev+1)
        field = .false.
      elseif (bfld.eq.'Q_STAR'.and.q_star)then
        write(2,*)'Q_STAR='
        write(2,'(10(E15.6,","))')((gradq_um(k,jt),jt=jt0,jt1),k=1,nlev+1)
        !   print *,'glr create gradq_um ',L_windrlx,((gradq_um(k,jt),jt=jt1-1,jt1-1),k=1,nlev+1)
        field = .false.
      endif

    endif

    if (bfld.eq.'UI'.and.ui)then
      write(2,*)'UI='
      write(2,'(10(E15.6,","))')(u_um(k,jt0),k=1,umlev)
      field = .false.
    elseif(bfld.eq.'VI'.and.vi)then
      write(2,*)'VI='
      write(2,'(10(E15.6,","))')(v_um(k,jt0),k=1,umlev)
      field = .false.

    elseif(bfld.eq.'WI'.and.wi)then
      write(2,*)'WI='
      write(2,'(10(E15.6,","))')(w_um(k,jt0),k=1,umlev)
      field = .false.

    elseif(bfld.eq.'THETA'.and.theta)then
      write(2,*)'THETA='
      write(2,'(10(E15.6,","))')(pt_um(k,jt0),k=1,umlev+1)
      write(6,'(10(E15.6,","))')(pt_um(k,jt0),k=1,umlev+1)
      field = .false.
    elseif(bfld.eq.'QI'.and.qi)then
      write(2,*)'QI='
      write(2,'(10(E15.6,","))')(q_um(k,jt0),k=1,umlev+1)
      field = .false.
    elseif(bfld.eq.'P_IN'.and.p_in)then
      write(2,*)'P_IN='
      write(2,'(10(E15.6,","))')(p_um(k,jt0),k=1,umlev)
      field = .false.

    elseif(bfld.eq.'LAT')then
      write(2,*)'LAT=     ',usrlat,','
      field = .false.
    elseif(bfld.eq.'LONG')then
      write(2,*)'LONG=     ',usrlon,','
      field = .false.
    elseif(bfld.eq.'YEAR_INIT')then
      write(2,*)'YEAR_INIT=     ',year,','
      field = .false.
    elseif(bfld.eq.'MONTH_INIT')then
      write(2,*)'MONTH_INIT=     ',month,','
      field = .false.
    elseif(bfld.eq.'DAY_INIT')then
      write(2,*)'DAY_INIT=     ',day,','
      field = .false.
    elseif(bfld.eq.'HOUR_INIT')then
      write(2,*)'HOUR_INIT=     ',shour,','
      field = .false.
    elseif(bfld.eq.'NMININ')then
      write(2,*)'NMININ=     ',int((nfor-1)*dt*60),','	! -301=-5hrs1min !glr
      field = .false.
      20       elseif(bfld.eq.'TSTAR_FORCING')then
      ch7 = '*288.0,'
      21        write(2,*)'TSTAR_FORCING=     '
      22	  write(2,'(I2,A7)')nfor,ch7
      23        field = .false.

    else

      if (field)then
        write(2,*)trim(base)
        field = .false.
      endif

    endif

    89	 continue

    if (diag)then
      write(2,*)trim(base)
      diag = .false.
    endif

  enddo

  9	continue

  print *,'NFOR = ',1+jt1-jt0
  !        if (L_windrlx) print *,'NFOR = ',1+jt1-jt0
  !        if (.not.L_windrlx) print *,'NFOR = ',jt1-jt0

  print *,''
  print *,'!-----------------------------------------------------------------------'
  print *,'! nb: the ichgf setting is sensitive to timestep and forcing frequency'
  print *,'!     and needs adjusting according to:'
  print *,'!'
  print *,'!     ichgf = (Time between forcing profiles) / (Timestep length)'
  print *,'!              [units = sec]'
  print *,'!-----------------------------------------------------------------------'
  print *,''
  print *,'================================'
  print *,'  genesis namelist.scm written'
  print *,'================================'
  print *,''
  print *,'Ensure you check through namelist settings before running SCM'
  print *,''



end subroutine create_namelist


