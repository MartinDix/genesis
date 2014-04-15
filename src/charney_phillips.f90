
	subroutine charney_phillips(nrec,nlev,levs,z,u,v,w,t,q,pt,msl,gradt,gradq,  &
     &	u_um,v_um,w_um,t_um,q_um,p_um,pt_um,gradt_um,gradq_um,debug)

!	subroutine charney_phillips(nrec,nlev,levs,z,u,ug,v,vg,w,t,q,pt,msl,gradt,gradq,  &
!     &	u_um,ug_um,v_um,vg_um,w_um,t_um,q_um,p_um,pt_um,gradt_um,gradq_um,debug)
	
!
!	Interpolates NWP fields onto UM Charney-Phillips grid
!	> adapted from Peter H's tapmscm2scum program
!								- vjb 7/5/09

	use global
	
	implicit none
	
!------ Input variables ------
!

	integer,intent(in)			:: nlev,nrec
	real,dimension(nlev),intent(in)		:: levs
	real,dimension(nlev,nrec),intent(in)	:: z,u,v,t,q,pt,w
!	real,dimension(nlev,nrec),intent(in)	:: z,u,v,t,q,pt,ug,vg,w
	real,dimension(nrec),intent(in)		:: msl
	real,dimension(nlev,nrec),intent(in)	:: gradt,gradq
	logical,intent(in)			:: debug


!------ Working variables and arrays ------
!
	integer				:: i,j,k,rec
	integer				:: kk
	integer				:: lwr,upr
	real				:: fact
	real				:: dp
	real				:: cd
	real				:: usfc
	real				:: ustar
	real				:: wspd
	real				:: aa,ab

	real,dimension(:),allocatable	:: zzr
	real,dimension(:),allocatable	:: zzt
	real,dimension(nlev,nrec)	:: relh
	real,dimension(nrec)		:: n
	


!====== Extras ======!
!------ Base namelist variables ------
!
	integer				:: nzum
	integer				:: sdate,shour
	integer				:: edate,ehour
	real				:: z_terrain_asl
	

!------ Vertlevs namelist variables ------
!
	integer				:: first_constant_r_rho_level
	real				:: z_top_of_model
	real,dimension(umlev+1)		:: eta_theta
	real,dimension(umlev)		:: eta_rho
	real,dimension(:),allocatable	:: zrho
	real,dimension(:),allocatable	:: ztheta	   !! UM level heights


!------ Output ------
!
	real,dimension(umlev+1,nrec),intent(out)	:: p_um
	real,dimension(umlev,nrec),intent(out)		:: u_um,v_um
!	real,dimension(umlev,nrec),intent(out)		:: u_um,v_um,ug_um,vg_um  ! (x)
	real,dimension(umlev+1,nrec),intent(out)	:: w_um,t_um,q_um,pt_um
	real,dimension(umlev+1,nrec),intent(out)	:: gradt_um,gradq_um

	
!------ Namelists ------
!
	namelist/base/z_terrain_asl,nzum
	namelist/vertlevs/z_top_of_model,first_constant_r_rho_level, &
     &                    eta_theta,eta_rho


!------ Read namelists ------
!
	read(*,base)
	read(*,vertlevs)

!	print *,'...OK'
!        print *,'glr char z_terrain_asl,nzum',z_terrain_asl,nzum
!        print *,'glr char z_top,frst ',z_top_of_model,first_constant_r_rho_level
!        print *,'glr char eta_theta ',eta_theta
!        print *,'glr char eta_rho ',eta_rho

!====================!     

	allocate(zrho(umlev+1))
	allocate(ztheta(umlev+1))
	allocate(zzt(umlev+1))
	allocate(zzr(umlev))
		
!====================!	
	
!------ Set UM level heights ------
!
	do k=1,umlev+1  
	  ztheta(k) = eta_theta(k) * z_top_of_model
	  zzt(k) = ztheta(k) + z_terrain_asl
	enddo
	
	do k=1,umlev
	  zrho(k) = eta_rho(k) * z_top_of_model
	  zzr(k) = zrho(k) + z_terrain_asl
	enddo
			

!------ Find pressure on model levels ------
!

	do rec=1,nrec

	  do k=1,umlev

	    if (zzr(k).lt.z(1,nrec))then	! test for UM level ht < lowest input data level
						! nb: usually ht of 1000hPa lev
	      dp = -1*(rho*gravity)*zzr(k)	! use hydrostatic eqn to get pressure nr sfc
	      p_um(k,rec) = msl(rec)+dp

	    elseif (zzr(k).gt.z(nlev,rec))then    ! if data > z_top

	      fact = (zzr(k) - z(nlev,rec)) / (maxz - z(nlev,rec))
	      p_um(k,rec) = levs(nlev)*(1.0-fact)+fact*minp   ! nb: minp = 0 !

	    else				! linearly interpolate to UM levels

	      do i=1,nlev			! upward sweep to find interpolation bounds

	        if (zzr(k).gt.z(i,rec))then	! find vertical upper and lower bounds
	          continue
                else
	          goto 50
	        endif

	      enddo

50	      continue
	      lwr = i-1				! i = level upper and lower bounds index
	      upr = i

	      fact = (zzr(k)-z(lwr,rec))/(z(upr,rec)-z(lwr,rec))	! linear
	      p_um(k,rec) = levs(lwr)*(1.0-fact)+fact*levs(upr)		! interpolation

	    endif

	  enddo

	  p_um(umlev+1,rec) = 100.0		! upper bound for 'extra' UM level

	enddo

!        print *,'glr charney shape z,zzr,zzt=',shape(z),shape(zzr),shape(zzt)
!        print *,'glr charney.. z=',z
!        print *,'glr charney.. zzr=',zzr
!        print *,'glr charney.. zzt=',zzt
!	print *,'vjb charney.. msl = ',msl
!	print *,'vjb charney.. p_in = ',p_um

!	stop


!------ Derive ustar from lowest input data level (Venkatram (1980) ------
!
	cd = (vkman / log(zzr(k)/zrough))		! drag coefficient


!------ Interpolate variables onto UM rho levels ------
!

	do rec=1,nrec

	  do k=1,umlev

	    cd = (vkman / log(zzr(1)/zrough))		! drag coefficient

	    if (zzr(k).lt.z(1,rec))then				! log interpolation to sfc
	      							! (neutral/stable assumption)
	      wspd = ( u(1,rec)**2 + v(1,rec)**2 )**0.5

!	      usfc = (4.7 * z(1,rec)) / (0.4 * 1100.0)		! for the record, V(80) uses this..
	      usfc = 0.2					! ..but z(1) can be >> 10m
								! so i'm hard-wiring it for now
	      aa = (cd * wspd)
	      ab = (2 * usfc) / ((cd**0.5)*wspd)
	      ustar = aa * (0.5 + 0.5*(1 - ab**2)**0.5)

	      print *,'zzr(k),wspd,usfc,cd,aa,ab,ustar = ',zzr(k),wspd,usfc,cd,aa,ab,ustar

	      if (ustar.lt.0.08)then
	        ustar = 0.08				! low wind speed limit on ustar
	      endif

	      u_um(k,rec) = (ustar/vkman)*log(zzr(k)/zrough)	! neutral PBL assumption
	      v_um(k,rec) = (ustar/vkman)*log(zzr(k)/zrough)	! neutral PBL assumption

	    elseif (zzr(k).gt.z(nlev,rec))then				! interpolate to the top

	      fact = (zzr(k) - z(nlev,rec)) / (maxz - z(nlev,rec))
	      u_um(k,rec) = u(nlev,rec)*(1-fact) + fact*zero
	      v_um(k,rec) = v(nlev,rec)*(1-fact) + fact*zero

	    else

	      do i=1,nlev

	        if (zzr(k).gt.z(i,rec))then
	          continue
	        else
	          goto 60
	        endif

	      enddo

60	      continue
	      lwr = i-1
	      upr = i

	      fact = (zzr(k) - z(lwr,rec)) / (z(upr,rec) - z(lwr,rec))
	      u_um(k,rec) = u(lwr,rec)*(1-fact) + (fact*u(upr,rec))
	      v_um(k,rec) = v(lwr,rec)*(1-fact) + (fact*v(upr,rec))

	    endif

	  enddo

	enddo


!	print *,'vjb charney.. u_um = ',u_um
!	print *,'vjb charney.. v_um = ',v_um



!------ Derive relh from spec hum ------
!

	call q2relh(nlev,nrec,levs,t,q,relh)
	print *,'relh = ',relh


!------ Derive cloud amount from relh profile (Slingo, 1980) ------
!
	call slingo(nlev,nrec,relh,levs,n)

	print *,'n = ',n

	stop


!------ Interpolate variables onto UM theta levels ------
!

	  do k = 1,nzum+1 !glr

	    t_um(k,rec) = 273.0			! ...default
	    q_um(k,rec) = 0.0			!      values
	    pt_um(k,rec) = 400.0		!         if nothing returned
	    gradt_um(k,rec) = 0.0
	    gradq_um(k,rec) = 0.0
!glr test making default the lowest obs values
	    w_um(k,rec) = w(1,rec)			! ...default
	    t_um(k,rec) = t(1,rec)			! ...default
	    q_um(k,rec) = q(1,rec)			!      values
	    pt_um(k,rec) = pt(1,rec)		!         if nothing returned
	    gradt_um(k,rec) = gradt(1,rec)
	    gradq_um(k,rec) = gradq(1,rec)

	    
	    if (zzt(k).le.z(nlev,rec))then	! test UM height < top of data
	   					!	if true...continue
	      do kk=2,nlev
	        
	        if (zzt(k).ge.z(kk-1,rec).and.zzt(k).le.z(kk,rec))then
! find fraction
! 		print *,'z-1,zzt,z = ',k,z(kk-1,rec),zzt(k),z(kk,rec)
		  
		  fact = (zzt(k) - z(kk-1,rec)) / (z(kk,rec) - z(kk-1,rec))
		  w_um(k,rec)  =  w(kk-1,rec) * (1.0-fact) + fact * w(kk,rec)
		  t_um(k,rec)  =  t(kk-1,rec) * (1.0-fact) + fact * t(kk,rec)
		  q_um(k,rec)  =  q(kk-1,rec) * (1.0-fact) + fact * q(kk,rec)
		  pt_um(k,rec) = pt(kk-1,rec) * (1.0-fact) + fact *pt(kk,rec)
		  gradt_um(k,rec) = gradt(kk-1,rec)*(1.0-fact)+fact*gradt(kk,rec)
		  gradq_um(k,rec) = gradq(kk-1,rec)*(1.0-fact)+fact*gradq(kk,rec)
		  
		  goto 200
		  
		endif
		
	      enddo
	     
	     else
	      
	      if (rec.eq.1) print *,'Theta lev (',k,') no data, default values set'
	    endif
	    
200	    continue
	    
!	    if (debug.and.rec.eq.1) print *,'w_um(k,1) = ',fact,w(kk-1,1),w_um(k,1),w(kk,1)
!	    if (debug.and.rec.eq.1) print *,'gradt_um(k,1) = ',gradt_um(kk,1)

	  
	   enddo
	  

	return

	end subroutine charney_phillips
