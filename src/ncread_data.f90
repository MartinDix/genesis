subroutine ncread_data(ncid,nvar,nlon,nlat,nlev,nrec,  &
    &                          lon,lat,lev,rec,dat)

    !	fills arrays according to NetCDF dimensions
    !	- vjb 20/8/09

    use global

    !        use netcdf
    implicit none
    include 'netcdf.inc'

    integer						:: i,j,k,l
    integer						:: nclen
    integer,intent(in)				:: ncid
    integer,intent(in)				:: nvar

    integer,intent(in)				:: nlon
    integer,intent(in)				:: nlat
    integer,intent(in)				:: nlev
    integer,intent(in)				:: nrec

    real,dimension(nlon),intent(inout)		:: lon
    real,dimension(nlat),intent(inout)		:: lat
    real,dimension(nlev),intent(inout)		:: lev
    real,dimension(nrec),intent(inout)		:: rec

    real,dimension(nlon,nlat,nlev,nrec),intent(inout)	:: dat

    character*20					:: var(num)
    character*20					:: name

    integer						:: lat_varid
    integer						:: lon_varid
    integer						:: lvl_varid
    integer						:: rec_varid
    integer						:: dat_varid

    logical                                         :: dd


    dd = .false.

    print *,'ncread_data: shape(dat) = ',shape(dat)


    if (nvar.eq.4)then

        dd = .true.
        if (dd) print *,'ncread_data: reduced dimension flag raised'
        do i=1,nvar
            call check(nf_inq_varname(ncid,i,name))
            var(i) = name
            print *,'var = ',var(i)

            if (i.eq.1)then
                call check(nf_inq_varid(ncid,var(i),lon_varid))
                call check(nf_get_var_real(ncid,lon_varid,lon))

            elseif (i.eq.2)then
                call check(nf_inq_varid(ncid,var(i),lat_varid))
                call check(nf_get_var_real(ncid,lat_varid,lat))

            elseif (i.eq.3)then
                call check(nf_inq_varid(ncid,var(i),rec_varid))
                call check(nf_get_var_real(ncid,rec_varid,rec))

            elseif (i.eq.4)then
                call check(nf_inq_varid(ncid,var(i),dat_varid))
                call check(nf_get_var_real(ncid,dat_varid,dat))

            endif

            lev = 1000.0

            print *,'read_data iteration = ',i

        enddo

    else

        do i=1,nvar
            call check(nf_inq_varname(ncid,i,name))
            var(i) = name
            print *,'var = ',var(i)

            if (i.eq.1)then
                call check(nf_inq_varid(ncid,var(i),lon_varid))
                call check(nf_get_var_real(ncid,lon_varid,lon))

            elseif (i.eq.2)then
                call check(nf_inq_varid(ncid,var(i),lat_varid))
                call check(nf_get_var_real(ncid,lat_varid,lat))

            elseif (i.eq.3)then
                call check(nf_inq_varid(ncid,var(i),lvl_varid))
                call check(nf_get_var_real(ncid,lvl_varid,lev))

            elseif (i.eq.4)then
                call check(nf_inq_varid(ncid,var(i),rec_varid))
                call check(nf_get_var_real(ncid,rec_varid,rec))

            elseif (i.eq.5)then
                call check(nf_inq_varid(ncid,var(i),dat_varid))
                call check(nf_get_var_real(ncid,dat_varid,dat))

            endif

            print *,'read_data iteration = ',i

        enddo

    endif

    call check(nf_close(ncid))


    return

  !------ Subroutines ------
  !  NetCDF function

contains

    subroutine check(status)
        integer, intent ( in) :: status

        if(status /= nf_noerr) then
            print *, trim(nf_strerror(status))
            stop "Stopped"
        end if
    end subroutine check


end subroutine ncread_data


