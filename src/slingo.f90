
	subroutine slingo(nlev,nrec,relh,levs,n)


	use global

	implicit none


!------ in/out parameters ------
!

	real,dimension(nlev,nrec),intent(in)		:: relh
	real,dimension(nlev),intent(in)			:: levs
	integer,intent(in)				:: nlev,nrec

!------ counters -------
!
	integer						:: k,rec
	integer						:: lyrs1,lyrs2,lyrs3

!------ working variables ------
!
	real						:: nl,nm,nh

!------ output -------
!
	real,dimension(nrec),intent(out)		:: n


!====== kick off =========
!


!------ find no. levels per layer for averaging ------
!
	lyrs1 = 0
	lyrs2 = 0
	lyrs3 = 0

	do k=1,nlev

	  if (levs(k).ge.80000)then
	    lyrs1 = lyrs1 + 1
	  elseif (levs(k).lt.80000.and.levs(k).ge.37000)then
	    lyrs2 = lyrs2 + 1
	  elseif (levs(k).lt.37000.and.levs(k).ge.12500)then
	    lyrs3 = lyrs3 + 1
	  endif

	enddo

	do rec=1,nrec

	  do k=1,nlev

	    if (levs(k).ge.80000)then

	      if (relh(k,rec).ge.80.0)then
	        nl = ((relh(k,rec)-80)**2) / 400
	      else
	        nl = 0.0
	      endif

	    elseif (levs(k).lt.80000.and.levs(k).ge.37000)then

	      if (relh(k,rec).ge.65.0)then
	        nm = ((relh(k,rec)-65)**2) / 1225
	      else
	        nm = 0.0
	      endif

	    elseif (levs(k).lt.37000.and.levs(k).ge.12500)then

	      if (relh(k,rec).ge.80.0)then
	        nh = ((relh(k,rec)-80)**2) / 400
	      else
	        nh = 0.0
	      endif

	    endif

	  enddo

	n(rec) = nl + nm + nh
	
	if (n(rec).ge.1)then
	  n(rec) = 0.0
	endif

	enddo

	print *,'slingo: layers1,2,3 = ',lyrs1,lyrs2,lyrs3

	end subroutine slingo
