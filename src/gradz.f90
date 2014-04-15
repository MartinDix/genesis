
        subroutine gradz(nrec,nlev,levs,temp,spechum,omega,gradtz,gradqz)

!       Derives vertical gradients of temperature and moisture on pressure
!       half-levels for inclusion in advective tendency calculation. 
!       - vjb 13/1/2011

        use global

        implicit none

!------ Input ------
!
        integer,intent(in)				:: nrec           
        integer,intent(in)				:: nlev

        real,dimension(nlev),intent(in)			:: levs

        real,dimension(nlev,nrec),intent(in)		:: temp
        real,dimension(nlev,nrec),intent(in)		:: spechum
        real,dimension(nlev,nrec),intent(in)		:: omega

!------ Working variables/arrays ------
!
        integer						:: i,j,k,rec

        real,dimension(nlev)				:: dp
        real,dimension(nlev)				:: dpm1
        real,dimension(nlev)				:: dpp1

        real,dimension(nlev,nrec)			:: dt
        real,dimension(nlev,nrec)			:: dq


!------ Output arrays ------
!
        real,dimension(nlev,nrec)			:: gradtz
        real,dimension(nlev,nrec)			:: gradqz


!------ kick off ------
!

        do rec=1,nrec
         do k=1,nlev

          if (k.eq.1.and.omega(k,rec).gt.0)then

           dp(k) = (levs(k+1) - levs(k)) / 2

           dt(k,rec) = (temp(k+1,rec) - temp(k,rec)) / 2
           dq(k,rec) = (spechum(k+1,rec) - spechum(k,rec)) / 2

           gradtz(k,rec) = dt(k,rec) / dp(k)
           gradqz(k,rec) = dq(k,rec) / dp(k)

          elseif (k.eq.1.and.omega(k,rec).le.0)then

!           dp(k) = 0

           dt(k,rec) = 0
           dq(k,rec) = 0

           gradtz(k,rec) = 0
           gradqz(k,rec) = 0

          elseif (k.eq.nlev.and.omega(k,rec).le.0)then

	   dp(k) = (levs(k-1) - levs(k)) / 2

	   dt(k,rec) = (temp(k-1,rec) - temp(k,rec)) / 2
           dq(k,rec) = (spechum(k-1,rec) - spechum(k,rec)) / 2

	   gradtz(k,rec) = dt(k,rec) / dp(k)
           gradqz(k,rec) = dq(k,rec) / dp(k)

          elseif (k.eq.nlev.and.omega(k,rec).gt.0)then

!	   dp(k) = 0

           dt(k,rec) = 0
           dq(k,rec) = 0

           gradtz(k,rec) = 0
           gradqz(k,rec) = 0

	  else

           dp(k) = ((levs(k+1) + levs(k))/2) - ((levs(k-1) + levs(k))/2)

           dt(k,rec) = ((temp(k+1,rec) + temp(k,rec))/2) - ((temp(k-1,rec) + temp(k,rec))/2)

           dq(k,rec) = ((spechum(k+1,rec) + spechum(k,rec))/2) - ((spechum(k-1,rec) + spechum(k,rec))/2)

           gradtz(k,rec) = dt(k,rec) / dp(k)
           gradqz(k,rec) = dq(k,rec) / dp(k)

          endif

         enddo

        enddo 

	print *,'gradtz,gradqz =',gradtz,gradqz

        return

        end subroutine gradz
