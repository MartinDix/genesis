subroutine q2relh(nlev,nrec,lev,temp,hum,relh)

  !  Takes input temperature (K) and specific humidity (kg/kg) and calculates corresponding
  !  relative humidity (%)
  !
  !  - vjb 18/4/2011

  implicit none

  !------ Declarations ------
  !
  integer,intent(in)  :: nlev
  integer,intent(in)  :: nrec

  real,dimension(nlev),intent(in)  :: lev    ! in Pascals
  real,dimension(nlev,nrec),intent(in)  :: temp   ! temperature
  real  :: temp1  ! temp T

  real,dimension(nlev,nrec),intent(in)  :: hum    ! q
  real,dimension(nlev,nrec),intent(out)  :: relh   ! relh
  real  :: es1    ! sat. vap. press.
  real  :: vap    ! vap. press.

  integer  :: k,m

  !------ Kick off ------
  !
  do m=1,nrec
    do k=1,nlev

      temp1 = temp(k,m)
      call establ(es1,temp1)
      vap = (hum(k,m)*lev(k))/0.622
      relh(k,m) = (vap/es1)*100

    enddo
  enddo

  return

end subroutine q2relh
