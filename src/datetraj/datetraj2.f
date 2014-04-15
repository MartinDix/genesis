      program datetraj2
c
c * Author: Kevin Keay  Date: Apr 12 2009
c
c	- mods for genesis application (vjb 20/5/09)
c
c * Notes:
c
c * 22/04/2009: Fixed a bug with tinc= -24 and date0 ending in 00
c
      character*12 cdate0
      integer tinc, ndays 
      character*80 optarg
c * maxt is 10 days @ 1 hr = 240
      parameter (maxt=240)  ! Max. no. of time periods in trajectory
      character*12 cdatet(maxt)
      character*4800 chout  ! maxt * 20 - for output on screen
c
      if(iargc().eq.0)then
	write(*,*)'Usage: datetraj date0 tinc ndays'
	write(*,*)'  date0: origin date (yyyymmddhh) e.g 1997060100'
        write(*,*)'  tinc: time increment in whole hours e.g. 6'
        write(*,*)'  Note: Use a -ve value of tinc for',
     * ' backward trajectories'
        write(*,*)'  ndays: trajectory length in whole days e.g. 4'
        write(*,*)'Examples: datetraj 1997060100  6 4'
        write(*,*)'          datetraj 2000020112 -6 5'
	write(*,*)'The list of dates is written to: datetraj.2.txt'
cv	write(*,*)'The backward flag (1= backward) and start/end dates',
cv     *	' are written to: datetraj.3.txt'
	stop
      elseif(iargc().ne.3)then
	write(*,*)'ERROR: 3 arguments required'
	stop
      else
	call getarg (1,optarg)
	read(optarg,'(A)')cdate0
	call getarg (2,optarg)
	read(optarg,*)tinc
	call getarg (3,optarg)
	read(optarg,*)ndays
c
        iback= 0
        if(tinc.lt.0)then
	  iback= 1  ! Backward trajectories
	endif
c
        xtinc= abs(real(tinc)) ! Make tinc > 0 aand real
c
        call datetraj (cdate0,iback,xtinc,ndays,maxt,nt,cdatet)
c
ckk	write(*,*)'nt=',nt
	write(chout,'(240A)')(cdatet(k),k=1,nt)
cv        write(*,'(A)')chout(1:ilen(chout))
cv        open (3,file='datetraj.3.txt')
cv	write(3,*)iback,' ',cdatet(1),cdatet(nt)
cv	close(3)
      endif
c
      end

c ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
c ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((

      subroutine datetraj (cdate,iback,hr,ndays,maxt,nt,cdatet)
c
c * Author: Kevin Keay  Date: Feb 6 2001
c
c * Modified:
c
c   (1) 9/3/2009: Converted to a subroutine with a few minor changes
c
      character*12 cdate
      integer iback,ndays,ndt
      character*12 cdatet(maxt)
c
      character*12 cdate0
      integer ih0, dayinc, tinc
      character*80 optarg
      character*1200 chout
c
      open (2,file='datetraj.2.txt')
c
        cdate0= cdate
	read(cdate(9:),*)ih0
	if(iback.eq.0)then
          tinc= int(hr)
	else
          tinc= int(-hr)
	endif
c
        iy4= 1       ! 4-digit year
c
        ntimes= ndays*24/iabs(tinc)
c
	read(cdate0,'(I4,2I2)')iy0,im0,id0
        call datevalu (id0,im0,iy0,jul0)
c
        nt= 1
	call valudate (jul0,id0,im0,iy0)
	if(iy4.eq.1.)then
cv	  write(cdatet(nt),'(I4,3I2.2)')iy0,im0,id0,ih0
          write(cdatet(nt),'(I4,2I2.2,2x,I2)')iy0,im0,id0,ih0
	else
	  write(cdatet(nt),'(I2.2,3I2.2)')iy0,im0,id0,ih0
	endif
cv	write(2,*)nt,' ',cdatet(nt)
        write(2,*)cdatet(nt)
c
        call datevalu (id0,im0,iy0,jul0)
	jul= jul0
        ih= ih0
c
	do while (nt.le.ntimes)
	  nt= nt +1
	  ih= ih +tinc
	  if(iback.eq.0)then
	    if(ih.ge.24.)then
	      ih= mod(ih,24)
	      jul= jul +1  ! Next day
	    endif
	  else
	    if(ih.le.0)then
	      ih= 24. -mod(abs(ih),24)
	      if(ih.eq.24)then
	        ih= 0
                if(tinc.eq.-24)then ! Case of tinc= -24  and hh= 00
	          jul= jul -1
                endif
	      else
	        jul= jul -1 ! Previous day
	      endif
	    endif
	  endif
	  call valudate (jul,id,im,iy)
	  if(iy4.eq.1.)then
cv	    write(cdatet(nt),'(I4,3I2.2)')iy,im,id,ih
	    write(cdatet(nt),'(I4,2I2.2,2x,I2)')iy,im,id,ih
	  else
	    write(cdatet(nt),'(I2.2,3I2.2)')iy,im,id,ih
	  endif
cv	  write(2,*)nt,' ',cdatet(nt)
          write(2,*)cdatet(nt)
	enddo
c
      close(2)
c
      return
      end

c ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
c ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((

      subroutine convdate (date,day,month,year)
c
      integer date,day,month,year
c
      year = date/10000
      month = date/100 - year*100
      day = date - year*10000 - month*100
c
      return
      end

c ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
c (((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((( 

      subroutine valudate (dval,d,m,y)
c
c     gives day, month, and year for given Julian day
c        written by David Hooke 
c
c * Source: Rural Water Corporation (1994)
c
      integer  dval, d, m, y
      integer  dv, t0
c
      t0 = 4
      y = t0*400 + 1
      dv = dval
      y = y + (dv/146097)*400
      dv = mod (dv,146097)
      y = y + (dv/36524)*100
      if (dv .eq. 146096) go to 10
      dv = mod (dv,36524)
      y = y + (dv/1461)*4
      dv = mod (dv,1461)
      y = y + dv/365
      if (dv .eq. 1460) go to 10
      dv = mod (dv,365)
      m    = 1
c
      do while (dv .ge. monthlen(m,y))
           dv = dv - monthlen (m,y)
           m = m + 1
      enddo
c
      d = dv + 1
      go to 90
  10  continue
      y = y - 1
      m = 12
      d = 31
  90  continue
c
      return
      end

c ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
c ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((

      subroutine datevalu (d,m,y,dval)
c
c     gives julian day for day, month, and year
c       written by David Hooke
c
c * Source: Rural Water Corporation (1994)
c
      integer  d, m, y, dval, dy, t0
      t0 = 4
      dy = 0
c
c     Compute day of year
c
      do 1000 i = 1,m-1
           dy = dy + monthlen(i,y)     ! days at end of month (m-1)
1000     continue
cc      enddo
      dy = dy + d - 1                  ! days at date (with 1st Jan zero).
c
c     Compute days since start of accumulation period .....................
c       (Day zero is 01/01/1601, negative before then).
c
     0 dval = ((y-1)/100-(t0*4))*36524
     1        + (y-(t0*400+1))/400
     2        + mod ((y-1),100)*365
     3        + mod ((y-1),100)/4
     4        + dy
c
      if (y.le.t0*400) dval = dval - 1
c
      return
      end

c ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
c ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((

      integer function monthlen(imon,iyr)
c
c * Source: Rural Water Corporation (1994)
c
      dimension mday(12)
      data mday/31,28,31,30,31,30,31,31,30,31,30,31/
      monthlen = mday(imon)
      if(imon.eq.2) then
         if(mod(iyr,4).eq.0) then
            monthlen = 29
            if(mod(iyr,100).eq.0) then
               if(mod(iyr,400).ne.0) then 
                  monthlen = 28
               endif
            endif
         endif
      endif
c
      return
      end

c ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
c ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((

      subroutine getyday (day,month,year,yday)
c
c * Module: GETYDAY 
c
c * Purpose: Returns the year day (1-366) for a given date (DAY,MONTH,YEAR).
c
c * Author: Kevin Keay     Date: 9/6/94
c
      implicit none
      integer*4 day,month,year,yday
      integer*4 monlen(12)
      integer*4 i
      data monlen /31,28,31,30,31,30,31,31,30,31,30,31/
      if((4*(year/4).eq.year.and.100*(year/100).ne.year)
     &  .or.(400*(year/400).eq.year))then
	monlen(2)= 29
      else
	monlen(2)= 28
      endif
      if(month.eq.1)then
	yday= day
      else
	yday= 0
	do i=1,month-1
	  yday= yday +monlen(i)
	end do
	yday= yday +day
      endif
c
      return
      end

c ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
c ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((

      integer function ilen(string)
c
      character*1 c
      character*(*) string
c
      do i=len(string),1,-1
        c= string(i:i)
        if (c.ne.' ') goto 10
      enddo
c * String is wholly blank
      ilen= -1
      return
c
   10 ilen= i
c
      return
      end

c ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
c ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
