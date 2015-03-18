subroutine ncread_dim(ncid,nvar,nlon,nlat,nlev,nrec,  &
    &	                     debug,infile)

    !	First read of NetCDF to get array dimensions
    !	- vjb 23/4/09

    use global

    !        use netcdf
    implicit none
    include 'netcdf.inc'


    !------ Passed variables
    !
    logical,intent(in)			:: debug
    character*25,intent(in)			:: infile

    integer,intent(out)			:: nlon		!
    integer,intent(out)			:: nlat		! passed
    integer,intent(out)			:: nlev		! arrays 2
    integer,intent(out)			:: nrec		!

    integer					:: i,j,k,l,m	! counters

    !------ NetCDF I/O
    !
    integer,intent(out)			:: ncid
    integer					:: ndim
    integer,intent(out)			:: nvar
    integer          			:: nvari
    integer					:: nattr,unlimdimid

    integer					:: lvl_dimid
    integer					:: lat_dimid
    integer					:: lon_dimid
    integer					:: rec_dimid

    integer					:: lat_varid
    integer					:: lon_varid
    integer					:: lvl_varid
    integer					:: rec_varid
    integer					:: dat_varid

    character (len = nf_max_name)		:: name
    integer					:: nclen
    character*20				:: var(num)
    integer					:: length(num)

    logical                                 :: dd


    !------ Open file and fill arrays ------
    !

    dd = .false.

    if (debug) print *,'Reading: ',infile

    call check(nf_open(infile,NF_NOWRITE,ncid))

    call check(nf_inq(ncid,ndim,nvar,nattr,unlimdimid))

    print *,'nvar',nvar

    if (nvar.eq.4)then
        do i=1,nvar
            dd = .true.
            if (dd) print *,'ncread_dim:  reduced dimension flag raised'
            call check(nf_inq_varname(ncid,i,name))
            var(i) = name
            print *,'var = ',var(i)
            if (i.eq.1) call check(nf_inq_varid(ncid,var(i),lon_varid))
            if (i.eq.2) call check(nf_inq_varid(ncid,var(i),lat_varid))
            if (i.eq.3) call check(nf_inq_varid(ncid,var(i),rec_varid))  ! problem with ncid
            if (i.eq.4) call check(nf_inq_varid(ncid,var(i),dat_varid))
        enddo

      !          var(4) = var(3)
      !          var(5) = var(4)
      !          var(3) = "lvl"

    elseif(nvar.eq.5)then

        do i=1,nvar
            call check(nf_inq_varname(ncid,i,name))
            var(i) = name
            print *,'var = ',var(i)
            if (i.eq.1) call check(nf_inq_varid(ncid,var(i),lon_varid))
            if (i.eq.2) call check(nf_inq_varid(ncid,var(i),lat_varid))
            if (i.eq.3) call check(nf_inq_varid(ncid,var(i),lvl_varid))
            if (i.eq.4) call check(nf_inq_varid(ncid,var(i),rec_varid))
            if (i.eq.5) call check(nf_inq_varid(ncid,var(i),dat_varid))
        enddo

    else

        stop'ERROR: nvar dimension outside acceptable threshold'

        endif

        if (debug) print*,'..1..'

        if (dd)then

            do i=1,nvar
                if (i.eq.1) call check(nf_inq_dim(ncid,lon_varid,  &
                    &                           var(i),length(i)))
                if (i.eq.2) call check(nf_inq_dim(ncid,lat_varid,  &
                    &                           var(i),length(i)))
                !         if (i.eq.3) length(i) = 1

                if (i.eq.3)then
                    call check(nf_inq_dim(ncid,rec_varid,var(i),length(i)))
                    var(i+1) = var(i)                    ! rearrange data dimensions
                    var(i) = "lvl"                       ! to standardise format
                    length(i+1) = length(i)
                    length(i) = 1
                endif

                if (debug) print *,'var,length = ',var(i),length(i)
            enddo

        else

            do i=1,nvar-1
                if (i.eq.1) call check(nf_inq_dim(ncid,lon_varid,  &
                    &                           var(i),length(i)))
                if (i.eq.2) call check(nf_inq_dim(ncid,lat_varid,  &
                    &                           var(i),length(i)))
                if (i.eq.3) call check(nf_inq_dim(ncid,lvl_varid,  &
                    &                           var(i),length(i)))
                if (i.eq.4) call check(nf_inq_dim(ncid,rec_varid,  &
                    &                           var(i),length(i)))
                if (debug) print *,'var,length = ',var(i),length(i)
            enddo

        endif

        if (debug) print*,'..2..'

        nlon = length(1)
        print *,'nlon = ',nlon
        nlat = length(2)
        print *,'nlat = ',nlat
        nlev = length(3)
        print *,'nlev = ',nlev
        nrec = length(4)
        print *,'nrec = ',nrec

        dd = .false.


        !------- return
        !
        print *,"**** Successfully read input file: ",infile

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

        !-------------------------------------------
        !  flip 1d array in vertical

        !	subroutine arrswapk(nx,ny,nz,nrec,infld)
        subroutine arrswapk(nz,infld)

            integer nx,ny,nz,nrec
            !	real, intent(inout)		:: infld(nx,ny,nz,nrec)
            real, intent(inout)		:: infld(nz)
            real				:: outfld(nz)

            integer i,j,k,l
            integer kinv
            real 				:: swap


            do k=1,nz
                kinv = nz+1 - k
                outfld(k) = infld(kinv)
            enddo
            do k=1,nz
                infld(k) = outfld(k)
            enddo

            return

        end subroutine arrswapk

        !----------------------------------------------
        !  flip 3d array in vertical - repeat for nrecs

        subroutine arrswap3dk(nx,ny,nz,nrec,infld)

            integer 			:: nx,ny,nz,nrec
            real, intent(inout)		:: infld(nx,ny,nz,nrec)
            real				:: outfld(nx,ny,nz,nrec)

            integer				:: i,j,k,l
            integer				:: kinv
            real 				:: swap

            do j=1,nrec
                do k=1,nz
                    kinv = nz+1 - k
                    do l=1,ny
                        do m=1,nx
                            outfld(m,l,k,j) = infld(m,l,kinv,j)
                        enddo
                    enddo
                enddo
            enddo

            do j=1,nrec
                do k=1,nz
                    do l=1,ny
                        do m=1,nx
                            infld(m,l,k,j) = outfld(m,l,k,j)
                        enddo
                    enddo
                enddo
            enddo

            return

        end subroutine arrswap3dk

      !------------------------------------------

    end subroutine ncread_dim


