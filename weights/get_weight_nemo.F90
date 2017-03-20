program interpol
use grid
use map
use netcdf

implicit none
   integer, parameter :: jpfld=8
   INTEGER , PARAMETER ::   jp_wndi = 1           ! index of 10m wind velocity (i-component) (m/s)    at T-point
   INTEGER , PARAMETER ::   jp_wndj = 2           ! index of 10m wind velocity (j-component) (m/s)    at T-point
   INTEGER , PARAMETER ::   jp_humi = 3           ! index of specific humidity               ( - )
   INTEGER , PARAMETER ::   jp_qsr  = 4           ! index of solar heat                      (W/m2)
   INTEGER , PARAMETER ::   jp_qlw  = 5           ! index of Long wave                       (W/m2)
   INTEGER , PARAMETER ::   jp_tair = 6           ! index of 10m air temperature             (Kelvin)
   INTEGER , PARAMETER ::   jp_prec = 7           ! index of total precipitation (rain+snow) (kg/m2/s)
   INTEGER , PARAMETER ::   jp_snow = 8           ! index of snow (solid precipitation)      (kg/m2/s)
   integer :: jp, nyear, nmonth, nday, ierror
   character(len=264) :: cn_dir
   TYPE(FLD_N), DIMENSION(jpfld) ::   slf_i       ! array of namelist informations on the fields to read
   TYPE(FLD_N) ::   sn_wndi, sn_wndj, sn_humi, sn_qsr       ! informations about the fields to be read
   TYPE(FLD_N) ::   sn_qlw , sn_tair, sn_prec, sn_snow      !   "                                 "
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf   ! structure of input fields (file informations, fields read)
   TYPE(WGT) :: refwgt                            ! weight and indices type storage

      NAMELIST/namsbc_core/ cn_dir ,            &
         &                  sn_wndi, sn_wndj, sn_humi  , sn_qsr ,           &
         &                  sn_qlw , sn_tair, sn_prec  , sn_snow
   !locals
   INTEGER :: varid          ! netcdf variable identifier
   INTEGER :: status, ncid
   INTEGER :: ji,jj,jk
   INTEGER :: dimx, dimy, dimz, dimx2, dimy2
   INTEGER :: ddims(2)
   INTEGER, allocatable, DIMENSION(:,:)          :: bufi
   DOUBLE PRECISION, allocatable, DIMENSION(:,:) :: bufd
   CHARACTER (len=15) :: varname
   !
      !!---------------------------------------------------------------------
! open ocean bathy file

   status=NF90_OPEN ( 'bathy_meter.nc', NF90_NOWRITE, ncid )
   status=NF90_INQ_VARID  (ncid, 'nav_lon', varid)
   status=NF90_INQUIRE_VARIABLE  (ncid, varid, dimids=ddims)

      !!---------------------------------------------------------------------
! find dimensions of NEMO grid and allocate

   status=NF90_INQUIRE_DIMENSION (ncid, ddims(1), len=jpi)
   status=NF90_INQUIRE_DIMENSION (ncid, ddims(2), len=jpj)
   allocate(glamt(jpi,jpj),gphit(jpi,jpj))
   allocate(bufi(jpi,jpj),bufd(jpi,jpj))

      !!---------------------------------------------------------------------
! read the lon/lat of ocean grid and close

   status=NF90_GET_VAR(ncid, varid, glamt)
   status=NF90_INQ_VARID  (ncid, 'nav_lat', varid)
   status=NF90_GET_VAR(ncid, varid, gphit)
   status=NF90_CLOSE( ncid )

      !!---------------------------------------------------------------------

         ! set file information (default values)
         cn_dir = './'       ! directory in which the model is executed
         !
         ! (NB: frequency positive => hours, negative => months)
         !            !    file     !  variable  !  clim   ! 'yearly' or !
         !            !    name     !   name     !  (T/F)  !  'monthly'  !
         sn_wndi = FLD_N( 'uwnd10m' ,  'u_10'    , .false. ,   'yearly' )
         sn_wndj = FLD_N( 'vwnd10m' ,  'v_10'    , .false. ,   'yearly' )
         sn_qsr  = FLD_N( 'qsw'     ,  'qsw'     , .false. ,   'yearly' )
         sn_qlw  = FLD_N( 'qlw'     ,  'qlw'     , .false. ,   'yearly' )
         sn_tair = FLD_N( 'tair10m' ,  't_10'    , .false. ,   'yearly' )
         sn_humi = FLD_N( 'humi10m' ,  'q_10'    , .false. ,   'yearly' )
         sn_prec = FLD_N( 'precip'  ,  'precip'  , .false. ,   'yearly' )
         sn_snow = FLD_N( 'snow'    ,  'snow'    , .false. ,   'yearly' )
         !
         OPEN( 1,file='namelist',status='old',action='read')
         READ( 1, namsbc_core )
         CLOSE(1)

         !
         !
         ! set sf structure
         ALLOCATE( sf(jpfld), STAT=ierror )
         DO jp= 1, jpfld
            ALLOCATE( sf(jp)%fnow(jpi,jpj) )
         END DO
         !
         ! store namelist information in an array
         slf_i(jp_wndi) = sn_wndi   ;   slf_i(jp_wndj) = sn_wndj
         slf_i(jp_qsr ) = sn_qsr    ;   slf_i(jp_qlw ) = sn_qlw
         slf_i(jp_tair) = sn_tair   ;   slf_i(jp_humi) = sn_humi
         slf_i(jp_prec) = sn_prec   ;   slf_i(jp_snow) = sn_snow
         ! fill sf with slf_i and control print
         CALL fld_fill( sf, slf_i, cn_dir, 'sbc_blk_core', 'flux formulattion for ocean surface boundary condition', 'namsbc_core' )

! set starting date
    nyear=2005
    nmonth=4
    nday=15

! find actual file name and initialize atmo grid
    CALL fld_clname( sf(1), nyear, nmonth, nday )
    CALL get_atmo_grid(sf(1))
    refwgt%numwgt = 4 ! only consider bilinear case
    ALLOCATE(refwgt%data_jpij(refwgt%numwgt,jpi,jpj))
    ALLOCATE(refwgt%data_wgt (refwgt%numwgt,jpi,jpj))
    CALL get_weight(refwgt)

!----------------------------------------------------------
!create output file
!----------------------------------------------------------

      status=NF90_CREATE( 'met_gem_weight.nc', NF90_CLOBBER, ncid )

      ! define dimensions
      refwgt%ddims(1)=nxgr
      refwgt%ddims(2)=nygr
      status=NF90_DEF_DIM( ncid, 'x', jpi, dimx )
      status=NF90_DEF_DIM( ncid, 'y', jpj, dimy )
      status=NF90_DEF_DIM( ncid, 'lon', refwgt%ddims(1), dimx2 )
      status=NF90_DEF_DIM( ncid, 'lat', refwgt%ddims(2), dimy2 )
      status=NF90_DEF_DIM( ncid, 'numwgt', refwgt%numwgt, dimz )

! create the lon/lat variable of the original atmo file
   IF ( grid_type == 1 ) THEN
      status=NF90_DEF_VAR( ncid, 'lon', NF90_DOUBLE, dimx2, varid )
      status=NF90_DEF_VAR( ncid, 'lat', NF90_DOUBLE, dimy2, varid )
   ELSEIF ( grid_type == 2 ) THEN
      ddims(1)=dimx2
      ddims(2)=dimy2
      status=NF90_DEF_VAR( ncid, 'nav_lon', NF90_DOUBLE, ddims, varid )
      status=NF90_DEF_VAR( ncid, 'nav_lat', NF90_DOUBLE, ddims, varid )
   ENDIF

! create the indices/weight pair for each corner
      ddims(1)=dimx
      ddims(2)=dimy
      DO jk=1,refwgt%numwgt
         varname = ' '
         WRITE(varname,'(a3,i2.2)') 'src',jk
         status=NF90_DEF_VAR( ncid, TRIM(varname), NF90_INT, ddims, varid )
         varname = ' '
         WRITE(varname,'(a3,i2.2)') 'wgt',jk
         status=NF90_DEF_VAR( ncid, TRIM(varname), NF90_DOUBLE, ddims, varid )
      ENDDO
      status=NF90_ENDDEF( ncid )

!----------------------------------------------------------
! fill in the values
!----------------------------------------------------------

! fill first the original grid axes
   IF ( grid_type == 1 ) THEN
      status=NF90_INQ_VARID(ncid, 'lon', varid)
      status=NF90_PUT_VAR(ncid, varid, x_axis(1:nxgr))
      status=NF90_INQ_VARID(ncid, 'lat', varid)
      status=NF90_PUT_VAR(ncid, varid, y_axis(1:nygr))
   ELSEIF ( grid_type == 2 ) THEN
      status=NF90_INQ_VARID(ncid, 'nav_lon', varid)
      status=NF90_PUT_VAR(ncid, varid, x_axis_2d)
      status=NF90_INQ_VARID(ncid, 'nav_lat', varid)
      status=NF90_PUT_VAR(ncid, varid, y_axis_2d)
   ENDIF

! fill now the indices/weight pairs
      DO jk=1,refwgt%numwgt
         varname = ' '
         WRITE(varname,'(a3,i2.2)') 'src',jk
         write(*,*) 'writing variable : ',TRIM(varname)
         status=NF90_INQ_VARID(ncid, TRIM(varname), varid)
         bufi(:,:) = refwgt%data_jpij(jk,:,:)
         status=NF90_PUT_VAR(ncid, varid, bufi)
         write(*,*) 'status put',status
         varname = ' '
         WRITE(varname,'(a3,i2.2)') 'wgt',jk
         write(*,*) 'writing variable : ',TRIM(varname)
         status=NF90_INQ_VARID(ncid, TRIM(varname), varid)
         bufd(:,:) = refwgt%data_wgt(jk,:,:)
         status=NF90_PUT_VAR(ncid, varid, bufd)
      ENDDO

status=NF90_CLOSE( ncid )

end

