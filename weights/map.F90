MODULE map
   !!======================================================================
   !!                       ***  MODULE  fldread  ***
   !! Ocean forcing:  read input field for surface boundary condition
   !!=====================================================================
   !! History :  9.0  !  06-06  (G. Madec) Original code
   !!                 !  05-08  (S. Alderson) Modified for Interpolation in memory
   !!                 !         from input grid to model grid
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   fld_read      : read input fields used for the computation of the
   !!                   surface boundary condition
   !!----------------------------------------------------------------------
   USE grid       ! ocean dynamics and tracers

   IMPLICIT NONE

   TYPE, PUBLIC ::   FLD_N      !: Namelist field informations
      CHARACTER(len = 256) ::   clname      ! generic name of the NetCDF flux file
      CHARACTER(len = 34)  ::   clvar       ! generic name of the variable in the NetCDF flux file
      LOGICAL              ::   ln_clim     ! climatology or not (T/F)
      CHARACTER(len = 7)   ::   cltype      ! type of data file 'daily', 'monthly' or yearly'
   END TYPE FLD_N

   TYPE, PUBLIC ::   FLD        !: Input field related variables
      CHARACTER(len = 256)            ::   clrootname   ! generic name of the NetCDF file
      CHARACTER(len = 256)            ::   clname       ! current name of the NetCDF file
      CHARACTER(len = 34)             ::   clvar        ! generic name of the variable in the NetCDF flux file
      LOGICAL                         ::   ln_clim      ! climatology or not (T/F)
      CHARACTER(len = 7)              ::   cltype       ! type of data file 'daily', 'monthly' or yearly'
      double precision , ALLOCATABLE, DIMENSION(:,:)   :: fnow  ! input fields interpolated to now time step
   END TYPE FLD

   TYPE         ::   WGT        !: Input weights related variables
      CHARACTER(len = 256)                    ::   wgtname      ! current name of the NetCDF weight file
      INTEGER , DIMENSION(2)                  ::   ddims        ! shape of input grid
      INTEGER                                 ::   numwgt       ! number of weights (4=bilinear, 16=bicubic)
      INTEGER, DIMENSION(:,:,:), POINTER      ::   data_jpij    ! array of source integers
      double precision, DIMENSION(:,:,:), POINTER  ::   data_wgt     ! array of weights on model grid
   END TYPE WGT

   !!----------------------------------------------------------------------
   !!        !  2007-06  (F. Dupont) space interpolation
   !!        !  2008-10  (F. Dupont) revised
   !!        !  unlike the weight based interpolation of above, this one is a simple
   !!        !  bilinear interpolation with no masking of land (beware!)
   !!----------------------------------------------------------------------
! netcdf variable ID:
   INTEGER, PUBLIC ::   idlon=0   !: longitude
   INTEGER, PUBLIC ::   idlat=0   !: latitude
   INTEGER, PUBLIC ::   idtime=0   !: time
! spatial and time axis   
   INTEGER, PUBLIC :: nxgr,nygr,nt ! dimension declared in the netcdf file
   INTEGER, PUBLIC ::  grid_type=0 !=1 1d axis, =2 2d axis
   INTEGER, PUBLIC ::  xtype=0 !=0 no need to get an extra row; =1 yes
   INTEGER, PUBLIC ::  ytype=0 !=0 no need to get an extra row; =1 yes
   double precision, ALLOCATABLE, DIMENSION(:), PUBLIC :: x_axis,y_axis,t_axis
   double precision, ALLOCATABLE, DIMENSION(:,:), PUBLIC :: x_axis_2d,y_axis_2d
   double precision :: x_axis_origin ! could be 0 (360) or -180 for a global data set
   LOGICAL, PRIVATE :: ln_atmo_grid_init_done=.FALSE. ! check if the grid has already been loaded
   double precision :: rsmall=1d-8

   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2006) 
   !! $Id: fldread.F90 1730 2009-11-16 14:34:19Z smasson $
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS


   SUBROUTINE fld_clname( sdjf, kyear, kmonth, kday, ldstop )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_clopn  ***
      !!
      !! ** Purpose :   update the file name and open the file
      !!
      !! ** Method  :   
      !!----------------------------------------------------------------------
      TYPE(FLD), INTENT(inout)           ::   sdjf     ! input field related variables
      INTEGER  , INTENT(in   )           ::   kyear    ! year value
      INTEGER  , INTENT(in   )           ::   kmonth   ! month value
      INTEGER  , INTENT(in   )           ::   kday     ! day value
      LOGICAL  , INTENT(in   ), OPTIONAL ::   ldstop   ! stop if open to read a non-existing file (default = .TRUE.)

      ! build the new filename if not climatological data
      IF( .NOT. sdjf%ln_clim ) THEN   ;   WRITE(sdjf%clname, '(a,"_y",i4.4)' ) TRIM( sdjf%clrootname ), kyear    ! add year
         IF( sdjf%cltype /= 'yearly' )    WRITE(sdjf%clname, '(a,"m" ,i2.2)' ) TRIM( sdjf%clname     ), kmonth   ! add month
         IF( sdjf%cltype == 'daily'  )    WRITE(sdjf%clname, '(a,"d" ,i2.2)' ) TRIM( sdjf%clname     ), kday     ! add day
      ELSE
      sdjf%clname = TRIM(sdjf%clrootname)
      ENDIF
      sdjf%clname = TRIM(sdjf%clname)//'.nc'
      !
   END SUBROUTINE fld_clname


   SUBROUTINE fld_fill( sdf, sdf_n, cdir, cdcaller, cdtitle, cdnam )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_fill  ***
      !!
      !! ** Purpose :   fill sdf with sdf_n and control print
      !!
      !! ** Method  :   
      !!----------------------------------------------------------------------
      TYPE(FLD)  , DIMENSION(:), INTENT(inout) ::   sdf        ! structure of input fields (file informations, fields read)
      TYPE(FLD_N), DIMENSION(:), INTENT(in   ) ::   sdf_n      ! array of namelist information structures
      CHARACTER(len=*)         , INTENT(in   ) ::   cdir       ! Root directory for location of flx files
      CHARACTER(len=*)         , INTENT(in   ) ::   cdcaller   ! 
      CHARACTER(len=*)         , INTENT(in   ) ::   cdtitle    ! 
      CHARACTER(len=*)         , INTENT(in   ) ::   cdnam      ! 
      !
      INTEGER  ::   jf       ! dummy indices
      !!---------------------------------------------------------------------

      DO jf = 1, SIZE(sdf)
         sdf(jf)%clrootname = TRIM( cdir )//TRIM( sdf_n(jf)%clname )
         sdf(jf)%clvar      = sdf_n(jf)%clvar
         sdf(jf)%ln_clim    = sdf_n(jf)%ln_clim
         sdf(jf)%cltype     = sdf_n(jf)%cltype
      END DO

      WRITE(*,*)
      WRITE(*,*) TRIM( cdcaller )//' : '//TRIM( cdtitle )
      WRITE(*,*) (/ ('~', jf = 1, LEN_TRIM( cdcaller ) ) /)
      WRITE(*,*) '          '//TRIM( cdnam )//' Namelist'
      WRITE(*,*) '          list of files'
      DO jf = 1, SIZE(sdf)
         WRITE(*,*) '               root filename: '  , TRIM( sdf(jf)%clrootname ),   &
            &                          ' variable name: '  , TRIM( sdf(jf)%clvar      ),   &
            &                          ' climatology: '    ,       sdf(jf)%ln_clim     ,   &
            &                          ' data type: '      ,       sdf(jf)%cltype
      END DO
      
   END SUBROUTINE fld_fill



   SUBROUTINE get_atmo_grid( sd )
      use netcdf
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  get_atmo_grid  ***
      !!
      !! ** Purpose : get the forcing atmospheric grid
      !!---------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE( FLD ),      INTENT(inout)         ::   sd            ! field with name of weights file
      !!
      INTEGER                                 ::   jn            ! dummy loop indices
      INTEGER                                 ::   id            ! temporary variable id
      CHARACTER (len=5)                       ::   aname
      INTEGER , DIMENSION(3)                  ::   ddims
      INTEGER , DIMENSION(jpi, jpj)           ::   data_src
      double precision, DIMENSION(jpi, jpj)           ::   data_tmp
      double precision, DIMENSION(:,:), ALLOCATABLE   ::   line2, lines  ! temporary array to read 2 lineumns
      CHARACTER (len=34)                      ::   lonvar, latvar
      LOGICAL                                 ::   cyclical
      double precision                                ::   resid, dlon   ! temporary array to read 2 lineumns
      INTEGER                                 ::   offset        ! temporary integer

      INTEGER, DIMENSION(3) ::   kdimsz    ! size of the dimensions

      CHARACTER(LEN=255)    ::   clname    ! the name of the file based on cdvar [[+clcpu]+clcpu]
      CHARACTER(LEN=128)    ::   varname   ! name of a variable in netcdf
      CHARACTER(LEN=100)    ::   clinfo    ! info character
      CHARACTER(LEN=4  )    ::   cyear     ! year in character
      CHARACTER(LEN=20)     ::   cdate   ! name of a variable in netcdf
      CHARACTER(LEN=2)      ::   cday    ! name of a variable in netcdf
      CHARACTER(LEN=3)      ::   cmonth2   ! name of a variable in netcdf
      LOGICAL               ::   llok      ! check the existence 
      INTEGER               ::   ncid, varid
      LOGICAL               ::   ll_fnd
      CHARACTER(LEN=100)    ::   time_origin,t_calendar ! calendar type declared in the netcdf file
      INTEGER               ::   ioipslid	   ! ioipsl file identifier
      INTEGER               ::   status
      INTEGER               ::   grid_dim(2) ! number of dimension for the forcing grid
      double precision              ::   xmin,xmax,dx,ymax
      INTEGER               ::   day,month,year,hour
      INTEGER               ::   nvars, ndims
      INTEGER               :: nxtmp, nytmp, nttmp
      double precision, ALLOCATABLE, DIMENSION(:) :: xtemp ! needed if xtype=1
      double precision, ALLOCATABLE, DIMENSION(:) :: ytemp ! needed if ytype=1
! netcdf variable ID:
      INTEGER               ::   idlon=0   !: longitude
      INTEGER               ::   idlat=0   !: latitude
      INTEGER               ::   idtime=0   !: time
      
      ! Initializations and control
      ! =============
      clinfo = '          get_atmo_grid ~~~  '

      !! open input data file (non-model grid)
      status=NF90_OPEN ( TRIM(sd%clname), NF90_NOWRITE, ncid )
      IF (status<0) THEN
       write(*,*) TRIM(clinfo)//' cannot find file '//TRIM(clname)
       stop
      ENDIF

! FD debug
write(*,*) 'reading : ',TRIM(sd%clname)
      grid_dim=0
      status=NF90_INQUIRE(ncid,ndims,nvars)
      DO varid = 1,nvars
        varname=''; ndims=0
        status=NF90_INQUIRE_VARIABLE(ncid, varid, name = varname, ndims = ndims)   ! number of dimensions
        IF (TRIM(varname) == "T" ) idtime = varid
        IF (TRIM(varname) == "TIME" ) idtime = varid
        IF (TRIM(varname) == "time" ) idtime = varid
        IF (varname(1:3) == "LON" .and. idlon==0 ) THEN; idlon = varid; grid_dim(1)=ndims; ENDIF
        IF (varname(1:3) == "lon" .and. idlon==0 ) THEN; idlon = varid; grid_dim(1)=ndims; ENDIF
        IF (TRIM(varname) == "nav_lon"           ) THEN; idlon = varid; grid_dim(1)=ndims; ENDIF
        IF (varname(1:3) == "LAT" .and. idlat==0 ) THEN; idlat = varid; grid_dim(2)=ndims; ENDIF
        IF (varname(1:3) == "lat" .and. idlat==0 ) THEN; idlat = varid; grid_dim(2)=ndims; ENDIF
        IF (TRIM(varname) == "nav_lat"           ) THEN; idlat = varid; grid_dim(2)=ndims; ENDIF
      ENDDO

      !! get dimensions
      ddims=0
      status=NF90_INQ_VARID         (ncid, sd%clvar, varid)
      status=NF90_INQUIRE_VARIABLE  (ncid, varid, dimids=kdimsz)
      IF (status.ne.0) THEN
       write(*,*) TRIM(clinfo)//' cannot find variable '//TRIM(sd%clvar)
       stop
      ENDIF
      ndims=0
      do id=1,3
         if (kdimsz(id)>0) ndims=ndims+1
      enddo
      IF (ndims.ne.3) THEN
         write(*,*) TRIM(clinfo)//' Dimensions should be 3 in variable '//TRIM(sd%clvar)
         stop
      ENDIF

      status=NF90_INQUIRE_DIMENSION (ncid, kdimsz(1), len=ddims(1))
      status=NF90_INQUIRE_DIMENSION (ncid, kdimsz(2), len=ddims(2))

      IF ( .not.ln_atmo_grid_init_done ) THEN
        nxgr = ddims(1); nygr = ddims(2)
        WRITE(*,*) 'atmospheric forcing netcdf grid dimensions: nx=',nxgr,', ny=',nygr
      ELSE
        IF (ddims(1).ne.nxgr) THEN
         write(*,*) TRIM(clinfo)//' Dimensions mismatch in x'
         stop
        ENDIF
        IF (ddims(2).ne.nygr) THEN
         write(*,*) TRIM(clinfo)//' Dimensions mismatch in y'
         stop
        ENDIF
      ENDIF

! return if the grid has already been loaded
      IF ( ln_atmo_grid_init_done ) THEN
         status=NF90_CLOSE( ncid ) ! close the file
         RETURN                    ! return to fldread
      ENDIF


      !! check for an east-west cyclic grid
      !! try to guess name of longitude variable

! FD had a return if finding the axis fails
      IF( idlon <= 0 ) THEN
       write(*,*) TRIM(clinfo)//' cannot find X axis variable '//TRIM(sd%clvar)
       stop
      ENDIF
      WRITE(*,*) TRIM(clinfo)//' found X axis varid:',idlon

! FD had a return if finding the axis fails
      IF( idlat <= 0 ) THEN
       write(*,*) TRIM(clinfo)//' cannot find Y axis variable '//TRIM(sd%clvar)
       stop
      ENDIF
      WRITE(*,*) TRIM(clinfo)//' found Y axis varid:',idlat

      offset = -1
      cyclical = .FALSE.
        !! get dimensions
      ddims=0
      status=NF90_INQUIRE_VARIABLE  (ncid, idlon, dimids=ddims)


      IF (grid_dim(1)==grid_dim(2)) THEN
        grid_type=grid_dim(1)
      ELSE
         write(*,*) TRIM(clinfo)//' incoherent axis type 1D/2D'
         stop
      ENDIF
      WRITE(*,*) 'grid_type',grid_type

      IF (ddims(2)==0) THEN  !! case where the axis are 1D


       !! found a longitude variable
       !! now going to assume that grid is regular so we can read a single row

       !! because input array is 2d, have to present iom with 2d array even though we only need 1d slice
       !! worse, we cant pass line2(:,1) to iom_get since this is treated as a 1d array which doesnt match input file
        ALLOCATE(xtemp(nxgr))
	  status=NF90_GET_VAR(ncid, idlon, xtemp)

       !! find largest grid spacing
        dlon = MAXVAL( xtemp(1:nxgr-1) )

        xmin=minval(xtemp)
        xmax=maxval(xtemp)
        x_axis_origin=xtemp(1)
        WRITE(*,'(A,3(1X,E13.6))') 'xmin/xmax/origin',xmin,xmax,x_axis_origin

        resid = ABS(ABS(xtemp(nxgr)-xtemp(1))-360.0)
        IF( resid < rsmall ) THEN
          !! end rows overlap in longitude
          offset = 0
          cyclical = .TRUE.
          xtype=0
          ALLOCATE(x_axis(nxgr))
          x_axis=xtemp
        ELSEIF( resid < 2.0*dlon ) THEN
          !! also call it cyclic if difference between end points is less than twice dlon from 360
          offset = 1
          xtype=1
          cyclical = .TRUE.
          ALLOCATE(x_axis(nxgr+1))
          x_axis(1:nxgr)=xtemp(1:nxgr)
!	  x_axis(nxgr+1)=mod(xtemp(1)+360.,360.)
          x_axis(nxgr+1)=xtemp(1)+360.
        ELSE ! local area domain
          xtype=2
          ALLOCATE(x_axis(nxgr))
          x_axis=xtemp
        ENDIF

        DEALLOCATE( xtemp )
        WRITE(*,*) 'xtype',xtype
        xmin=minval(x_axis)
        xmax=maxval(x_axis)
        x_axis_origin=x_axis(1)
        WRITE(*,'(A,3(1X,E13.6))') 'xmin/xmax/origin',xmin,xmax,x_axis_origin

! check if y_axis goes North enough to cover the oceanic domain
        ALLOCATE(ytemp(nygr))
        status=NF90_GET_VAR(ncid, idlat, ytemp)
        ymax=maxval(ytemp)
        WRITE(*,*) 'Max latitude in atmo. forcing file',ymax
        WRITE(*,*) 'Max latitude in model latitude    ',maxval(gphit)
        IF ( ymax-maxval(gphit) .lt. rsmall ) THEN
           ytype=1
	   ALLOCATE(y_axis(nygr+1))
           y_axis(1:nygr)=ytemp(1:nygr)
           y_axis(nygr+1)=90.
        ELSE
           ALLOCATE(y_axis(nygr))
           y_axis=ytemp
        ENDIF
        WRITE(*,*) 'ytype',ytype

        DEALLOCATE(ytemp)

      ELSE !! if the axes are 2D

        ALLOCATE(x_axis_2d(nxgr,nygr),y_axis_2d(nxgr,nygr))
        status=NF90_GET_VAR(ncid, idlon, x_axis_2d)
        status=NF90_GET_VAR(ncid, idlat, y_axis_2d)
! check if x_axis has a redundant 0/360
        xmin=minval(x_axis_2d)
        xmax=maxval(x_axis_2d)
        x_axis_origin=xmin
        WRITE(*,'(A,3(1X,E13.6))') 'xmin/xmax/origin',xmin,xmax,x_axis_origin

      ENDIF


      !! close it
      status=NF90_CLOSE( ncid )


      ln_atmo_grid_init_done = .true.

	 
   END SUBROUTINE get_atmo_grid
   

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


   SUBROUTINE get_data_cdf_2d(clname, varname, zdata, time_index)


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
   use netcdf
   IMPLICIT NONE
   !arguments
   CHARACTER(len=256), INTENT(IN) :: clname        ! netcdf file name
   CHARACTER(len=34), INTENT(IN) :: varname        ! netcdf variable name
   INTEGER, INTENT(IN)           :: time_index     ! netcdf file identifier
   double precision, DIMENSION(jpi,jpj), &
                    INTENT (OUT) :: zdata          ! buffer data
   !locals
   INTEGER                       :: varid          ! netcdf variable identifier
   INTEGER :: status, ncid
   INTEGER, DIMENSION(3) :: start,count
   double precision :: scale_factor, add_offset
   double precision, ALLOCATABLE, DIMENSION(:,:) :: buf, buf2
   double precision :: mask=-999.9
   !
   ALLOCATE(buf(nxgr,nygr))

   status=NF90_OPEN ( TRIM(clname), NF90_NOWRITE, ncid )
   status=NF90_INQ_VARID  (ncid, TRIM(varname), varid)
   count(1)=nxgr; count(2)=nygr; count(3)=1
   start(1)= 1; start(2)= 1; start(3)=time_index
! FD debug
!   write(*,*) 'debug',ncid, TRIM(varname), varid, start, count
   status=NF90_GET_VAR(ncid, varid, buf, start, count)
   status=NF90_CLOSE( ncid )


   IF ( grid_type == 2 ) THEN
     CALL map_interpolation_irreg(nxgr, nygr, x_axis_2d, y_axis_2d, buf, jpi, jpj, glamt, gphit, mask, zdata)
     DEALLOCATE(buf)
     RETURN
   ENDIF

   IF (xtype==1.and.ytype==0) THEN
      ALLOCATE(buf2(nxgr+1,nygr))
      buf2(1:nxgr,:)=buf(1:nxgr,:)
      buf2(nxgr+1,:)=buf(1,:)
      CALL map_interpolation_reg(nxgr+1, nygr, x_axis, y_axis, buf2, jpi, jpj, glamt, gphit, mask, zdata)
      DEALLOCATE(buf2)
!
   ELSE IF (xtype==1.and.ytype==1) THEN
      ALLOCATE(buf2(nxgr+1,nygr+1))
      buf2(1:nxgr,1:nygr)=buf(1:nxgr,1:nygr)
      buf2(nxgr+1,1:nygr)=buf(1,:nygr)
      buf2(:,nygr+1)=buf2(:,nygr)
      CALL map_interpolation_reg(nxgr+1, nygr+1, x_axis, y_axis, buf2, jpi, jpj, glamt, gphit, mask, zdata)
      DEALLOCATE(buf2)
!
   ELSE IF (xtype==0.and.ytype==1) THEN
      ALLOCATE(buf2(nxgr,nygr+1))
      buf2(:,1:nygr)=buf(:,1:nygr)
      buf2(:,nygr+1)=buf2(:,nygr)
      CALL map_interpolation_reg(nxgr, nygr+1, x_axis, y_axis, buf2, jpi, jpj, glamt, gphit, mask, zdata)
      DEALLOCATE(buf2)
!
   ELSE
      CALL map_interpolation_reg(nxgr, nygr, x_axis, y_axis, buf, jpi, jpj, glamt, gphit, mask, zdata)
   ENDIF

   DEALLOCATE(buf)

   END SUBROUTINE get_data_cdf_2d
   

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


subroutine map_interpolation_reg(nx, ny, xgr, ygr, zgr, &
nx2, ny2, xdest, ydest, mask, zdest)


!*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

implicit none

! arguments

integer nx, ny
integer nx2, ny2
double precision, dimension(nx, ny)   :: zgr
double precision, dimension(nx)       :: xgr
double precision, dimension(ny)       :: ygr
double precision, dimension(nx2, ny2) :: xdest, ydest, zdest
double precision                      :: mask

! local

double precision :: w,s,area(4),dx,dy,zz,x,y
integer  :: k,l,k0,l0
double precision :: ux,uy,un,vx,vy,vn
double precision :: xx(0:1),yy(0:1)
double precision :: value(4),zinf,zsup
integer  :: i,j,m,count,i0,j0,m0,i2,j2
integer  :: attempt,kprev,lprev,status
double precision :: a,b,d
double precision :: eps=1e-6
double precision :: grid_xmin,grid_xmax,grid_ymax,grid_ymin

grid_xmin=minval(xgr)
grid_xmax=maxval(xgr)
grid_ymin=minval(ygr)
grid_ymax=maxval(ygr)

do j2=1,ny2
 do i2=1,nx2
  x=xdest(i2,j2)
  y=ydest(i2,j2)
  
! for global dataset:
  if (x<x_axis_origin .and. xtype<=1 ) x=mod(x+360.,360.)

if (xtype==2) then
  if(x < grid_xmin) then
    write(*,*) 'Error in map_interpolation: xmin=',grid_xmin,' and x=',x
    return
  endif
  if (x > grid_xmax) then
    write(*,*) 'Error in map_interpolation: xmax=',grid_xmax,' and x=',x
    return
  endif
endif
  if(y < grid_ymin) then
    write(*,*) 'Error in map_interpolation: ymin=',grid_ymin,' and y=',y
    return
  endif
  if (y > grid_ymax) then
    write(*,*) 'Error in map_interpolation: ymax=',grid_ymax,' and y=',y
    return
  endif

  k0=nx/2
  l0=ny/2
  status=-1
  

  do attempt=1,10

    kprev=k0
    lprev=l0
    dx = xgr(k0+1)-xgr(k0)
    dy = ygr(l0+1)-ygr(l0)
    
    k0 = k0 + nint((x-xgr(k0))/dx-.5)
    if(k0<1) k0=1
    if(k0>=nx) k0=nx-1
    
    l0 = l0 + nint((y-ygr(l0))/dy-.5)
    if(l0<1) l0=1
    if(l0>=ny) l0=ny-1
    
!    write(*,*) attempt,k0,l0
    
    if((kprev==k0) .and. (lprev==l0)) exit

   enddo
   
! compute the a,b in first sector

    a=(x-xgr(k0))/dx
    b=(y-ygr(l0))/dy

! interpolate

  count=0

  xx(0)=1.-a
  xx(1)=a
    
  yy(0)=1.-b
  yy(1)=b

  do i=0,1
    do j=0,1
      zz=zgr(i+k0,j+l0)
      if (zz .ne. mask) then
	count = count + 1
	value(count)=zz
	area(count)=xx(i)*yy(j)
      endif
    enddo
  enddo

  w=0.0
  do m=1,count
    w=w+area(m)
  enddo

  zz=0.0;
  if (w > 0.0) then
    do m=1,count
      zz=zz+value(m)*area(m)
    enddo
    zz=zz/w
    zdest(i2,j2)=zz
  else
    zdest(i2,j2)=mask
  endif

 enddo
enddo

 return
end subroutine map_interpolation_reg



!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


subroutine map_interpolation_irreg(nx,ny, xgr, ygr, zgr, &
nx2, ny2, xdest, ydest, mask, zdest)


!*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

implicit none

! arguments

integer nx, ny
integer nx2, ny2
double precision, dimension(nx, ny) :: xgr,ygr,zgr
double precision, dimension(nx2, ny2) :: xdest, ydest, zdest
double precision :: mask

! local

double precision :: w,s,area(4),dx,dy,zz,x,y
integer k,l,k0,l0
double precision :: ux,uy,un,vx,vy,vn
double precision :: xx(0:1),yy(0:1)
double precision :: value(4),zinf,zsup
integer i,j,m,count,i0,j0,m0,i2,j2
integer attempt,kprev,lprev,status
integer redo_k,redo_l
double precision :: a,b,d
integer step_k,step_l
double precision :: eps=1e-6
double precision :: grid_xmin,grid_xmax,grid_ymax,grid_ymin

grid_xmin=minval(xgr)
grid_xmax=maxval(xgr)
grid_ymin=minval(ygr)
grid_ymax=maxval(ygr)


do j2=1,ny2
 do i2=1,nx2
  x=xdest(i2,j2)
  y=ydest(i2,j2)
  
  if (x<0) x=x+360
  if(x < grid_xmin) then
    write(*,*) 'Error in map_interpolation: xmin=',grid_xmin,' and x=',x
    return
  endif
  if (x > grid_xmax) then
    write(*,*) 'Error in map_interpolation: xmax=',grid_xmax,' and x=',x
    return
  endif
  if(y < grid_ymin) then
    write(*,*) 'Error in map_interpolation: ymin=',grid_ymin,' and y=',y
    return
  endif
  if (y > grid_ymax) then
    write(*,*) 'Error in map_interpolation: ymax=',grid_ymax,' and y=',y
    return
  endif

  k0=nx/2
  l0=ny/2
  status=-1
  
  do attempt=1,10

    kprev=k0
    lprev=l0

    ux = xgr(k0+1,l0)-xgr(k0,l0)
    uy = ygr(k0+1,l0)-ygr(k0,l0)

    vx = xgr(k0,l0+1)-xgr(k0,l0)
    vy = ygr(k0,l0+1)-ygr(k0,l0)

    dx=x-xgr(k0,l0)
    dy=y-ygr(k0,l0)

    d=ux*vy-uy*vx

    a=(dx*vy-dy*vx)/d
    b=(ux*dy-uy*dx)/d
    
    k0 = k0 + nint (a)
    l0 = l0 + nint (b)

    if(l0<1) l0=1
    if(l0>ny) l0=ny-1

    if(l0<1) l0=1
    if(l0>ny) l0=ny-1
    
!    write(*,*) attempt,k0,l0
    
    if((kprev==k0) .and. (lprev==l0)) exit

   enddo
   
! initial guess on the closest point

    ux = xgr(k0+1,l0)-xgr(k0,l0)
    uy = ygr(k0+1,l0)-ygr(k0,l0)

    vx = xgr(k0,l0+1)-xgr(k0,l0)
    vy = ygr(k0,l0+1)-ygr(k0,l0)

    dx=x-xgr(k0,l0)
    dy=y-ygr(k0,l0)

    d=ux*vy-uy*vx

    a=(dx*vy-dy*vx)/d
    b=(ux*dy-uy*dx)/d
    
!    write(*,*) a,b

    k0 = k0 + nint (a)
    l0 = l0 + nint (b)
    
!    write(*,*) k0,l0
  
  
! recompute the a,b in first sector

    ux = xgr(k0+1,l0)-xgr(k0,l0)
    uy = ygr(k0+1,l0)-ygr(k0,l0)

    vx = xgr(k0,l0+1)-xgr(k0,l0)
    vy = ygr(k0,l0+1)-ygr(k0,l0)

    dx=x-xgr(k0,l0)
    dy=y-ygr(k0,l0)

    d=ux*vy-uy*vx

    a=(dx*vy-dy*vx)/d
    b=(ux*dy-uy*dx)/d
  
!  write(*,*) '1 sec',a,b

  step_k=1
  step_l=1

  if (a>-eps .and. b>-eps) goto 1000

! recompute the a,b in second sector

    ux = xgr(k0-1,l0)-xgr(k0,l0)
    uy = ygr(k0-1,l0)-ygr(k0,l0)

    vx = xgr(k0,l0+1)-xgr(k0,l0)
    vy = ygr(k0,l0+1)-ygr(k0,l0)

    dx=x-xgr(k0,l0)
    dy=y-ygr(k0,l0)

    d=ux*vy-uy*vx

    a=(dx*vy-dy*vx)/d
    b=(ux*dy-uy*dx)/d
  
!  write(*,*) '2 sec',a,b

  step_k=-1
  step_l=1

  if (a>-eps .and. b>-eps) goto 1000

! recompute the a,b in third sector

    ux = xgr(k0-1,l0)-xgr(k0,l0)
    uy = ygr(k0-1,l0)-ygr(k0,l0)

    vx = xgr(k0,l0-1)-xgr(k0,l0)
    vy = ygr(k0,l0-1)-ygr(k0,l0)

    dx=x-xgr(k0,l0)
    dy=y-ygr(k0,l0)

    d=ux*vy-uy*vx

    a=(dx*vy-dy*vx)/d
    b=(ux*dy-uy*dx)/d
  
!  write(*,*) '3 sec',a,b

  step_k=-1
  step_l=-1

  if (a>-eps .and. b>-eps) goto 1000

! recompute the a,b in fourth sector

    ux = xgr(k0+1,l0)-xgr(k0,l0)
    uy = ygr(k0+1,l0)-ygr(k0,l0)

    vx = xgr(k0,l0-1)-xgr(k0,l0)
    vy = ygr(k0,l0-1)-ygr(k0,l0)

    dx=x-xgr(k0,l0)
    dy=y-ygr(k0,l0)

    d=ux*vy-uy*vx

    a=(dx*vy-dy*vx)/d
    b=(ux*dy-uy*dx)/d
  
!  write(*,*) '4 sec',a,b

  step_k=1
  step_l=-1

  if (a>-eps .and. b>-eps) then
   goto 1000
  else
   write(*,*) 'issues with interpolation in sectors'
   return
  endif

!-----------------------------------------------------------------------
! Interpolation (bilinear)

1000 continue

  count=0

  xx(0)=1.-a
  xx(1)=a
    
  yy(0)=1.-b
  yy(1)=b

  do i=0,1
    do j=0,1
      zz=zgr(step_k*i+k0,step_l*j+l0)
      if (zz .ne. mask) then
	count = count + 1
	value(count)=zz
	area(count)=xx(i)*yy(j)
      endif
    enddo
  enddo

  w=0.0
  do m=1,count
    w=w+area(m)
  enddo

  zz=0.0;
  if (w > 0.0) then
    do m=1,count
     zz=zz+value(m)*area(m)
    enddo
    zz=zz/w
    zdest(i2,j2)=zz
  else
    zdest(i2,j2)=mask
  endif

 enddo
enddo

 return
end subroutine map_interpolation_irreg


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


   SUBROUTINE get_weight(rwgt)


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
   use netcdf
   !!---------------------------------------------------------------------
   !!                   ***  SUBROUTINE  get_weight  ***
   !!
   !! ** Purpose : get the forcing atmospheric grid weighting
   !!---------------------------------------------------------------------
   IMPLICIT NONE
   !arguments
   TYPE( WGT ),INTENT(inout)  ::   rwgt          ! field with name of weights file
   !locals
   double precision :: mask=-999.9
   !
   IF ( grid_type == 2 ) THEN
     CALL weight_interpolation_irreg(nxgr, nygr, x_axis_2d, y_axis_2d, jpi, jpj, &
                glamt, gphit, mask, rwgt%data_jpij, rwgt%data_wgt)
     RETURN
   ENDIF

   IF (xtype==1.and.ytype==0) THEN
     CALL weight_interpolation_reg(nxgr+1, nygr, x_axis, y_axis, jpi, jpj, &
                glamt, gphit, mask, rwgt%data_jpij, rwgt%data_wgt)
!
   ELSE IF (xtype==1.and.ytype==1) THEN
     CALL weight_interpolation_reg(nxgr+1, nygr+1, x_axis, y_axis, jpi, jpj, &
                glamt, gphit, mask, rwgt%data_jpij, rwgt%data_wgt)
!
   ELSE IF (xtype==0.and.ytype==1) THEN
     CALL weight_interpolation_reg(nxgr, nygr+1, x_axis, y_axis, jpi, jpj, &
                glamt, gphit, mask, rwgt%data_jpij, rwgt%data_wgt)
!
   ELSE
     CALL weight_interpolation_reg(nxgr, nygr, x_axis, y_axis, jpi, jpj, &
                glamt, gphit, mask, rwgt%data_jpij, rwgt%data_wgt)
   ENDIF

   END SUBROUTINE get_weight


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


subroutine weight_interpolation_reg(nx, ny, xgr, ygr, &
nx2, ny2, xdest, ydest, mask, indexij, weight)


!*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

implicit none

! arguments

integer nx, ny
integer nx2, ny2
double precision, dimension(nx)       :: xgr
double precision, dimension(ny)       :: ygr
double precision, dimension(nx2, ny2) :: xdest, ydest
integer, dimension(4,nx2,ny2)         :: indexij
double precision, dimension(4,nx2,ny2):: weight
double precision                      :: mask

! local

double precision :: w,s,dx,dy,x,y
integer  :: k,l,k0,l0
double precision :: ux,uy,un,vx,vy,vn
double precision :: xx(0:1),yy(0:1)
integer  :: i,j,m,count,i0,j0,m0,i2,j2,i1,j1
integer  :: attempt,kprev,lprev,status
double precision :: a,b,d
double precision :: eps=1e-6
double precision :: grid_xmin,grid_xmax,grid_ymax,grid_ymin
integer :: nx0,ny0

! find true size of incoming forcing
nx0=nx; ny0=ny
IF (xtype==1) nx0=nx-1
IF (ytype==1) ny0=ny-1


grid_xmin=minval(xgr)
grid_xmax=maxval(xgr)
grid_ymin=minval(ygr)
grid_ymax=maxval(ygr)

do j2=1,ny2
 do i2=1,nx2
  x=xdest(i2,j2)
  y=ydest(i2,j2)
  
! for global dataset:
  if (x<x_axis_origin .and. xtype<=1 ) x=mod(x+360.d0,360.d0)
  if (x>360.d0) x=mod(x,360.d0)

if (xtype==2) then
  if(x < grid_xmin) then
    write(*,*) 'Error in map_interpolation: xmin=',grid_xmin,' and x=',x
    return
  endif
  if (x > grid_xmax) then
    write(*,*) 'Error in map_interpolation: xmax=',grid_xmax,' and x=',x
    return
  endif
endif
  if(y < grid_ymin) then
    write(*,*) 'Error in map_interpolation: ymin=',grid_ymin,' and y=',y
    return
  endif
  if (y > grid_ymax) then
    write(*,*) 'Error in map_interpolation: ymax=',grid_ymax,' and y=',y
    return
  endif

  k0=nx/2
  l0=ny/2
  status=-1
  

  do attempt=1,10

    kprev=k0
    lprev=l0
    dx = xgr(k0+1)-xgr(k0)
    dy = ygr(l0+1)-ygr(l0)
    
    k0 = k0 + nint((x-xgr(k0))/dx-.5)
    if(k0<1) k0=1
    if(k0>=nx) k0=nx-1
    
    l0 = l0 + nint((y-ygr(l0))/dy-.5)
    if(l0<1) l0=1
    if(l0>=ny) l0=ny-1
    
!    write(*,*) attempt,k0,l0
    
    if((kprev==k0) .and. (lprev==l0)) exit

   enddo
   
! compute the a,b in first sector

    a=(x-xgr(k0))/dx
    b=(y-ygr(l0))/dy

! interpolate

  xx(0)=1.-a
  xx(1)=a
    
  yy(0)=1.-b
  yy(1)=b

  count=0
  do i=0,1
    do j=0,1
      count = count + 1
      i1=i+k0; j1=j+l0
      IF (i1>nx0) i1=1
      IF (j1>ny0) j1=ny0
      indexij(count,i2,j2)=i1+nx0*(j1-1)
      weight(count,i2,j2)=xx(i)*yy(j)
! FD debug
if (weight(count,i2,j2).lt.0.d0) then
 write(*,*) 'weigth',i2,j2,k0,xgr(k0),xgr(k0+1),x
 stop
endif
    enddo
  enddo

  w=0.d0
  do m=1,count
    w=w+weight(m,i2,j2)
  enddo

  if (w > 0.d0) then
    do m=1,count
      weight(m,i2,j2)=weight(m,i2,j2)/w
    enddo
  else
    do m=1,count
      weight(m,i2,j2)=mask
    enddo
  endif

 enddo
enddo

 return
end subroutine weight_interpolation_reg



!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


subroutine weight_interpolation_irreg(nx,ny, xgr, ygr, &
nx2, ny2, xdest, ydest, mask, indexij, weight)


!*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

implicit none

! arguments

integer nx, ny
integer nx2, ny2
double precision, dimension(nx, ny) :: xgr,ygr
double precision, dimension(nx2, ny2) :: xdest, ydest
integer, dimension(4,nx2,ny2)         :: indexij
double precision, dimension(4,nx2,ny2):: weight
double precision :: mask

! local

double precision :: w,s,dx,dy,x,y
integer k,l,k0,l0
double precision :: ux,uy,un,vx,vy,vn
double precision :: xx(0:1),yy(0:1)
integer i,j,m,count,i0,j0,m0,i2,j2,i1,j1
integer attempt,kprev,lprev,status
integer redo_k,redo_l
double precision :: a,b,d
integer step_k,step_l
double precision :: eps=1e-6
double precision :: grid_xmin,grid_xmax,grid_ymax,grid_ymin

grid_xmin=minval(xgr)
grid_xmax=maxval(xgr)
grid_ymin=minval(ygr)
grid_ymax=maxval(ygr)


do j2=1,ny2
 do i2=1,nx2
  x=xdest(i2,j2)
  y=ydest(i2,j2)
  
  if (x<0) x=x+360
  if(x < grid_xmin) then
    write(*,*) 'Error in map_interpolation: xmin=',grid_xmin,' and x=',x
    return
  endif
  if (x > grid_xmax) then
    write(*,*) 'Error in map_interpolation: xmax=',grid_xmax,' and x=',x
    return
  endif
  if(y < grid_ymin) then
    write(*,*) 'Error in map_interpolation: ymin=',grid_ymin,' and y=',y
    return
  endif
  if (y > grid_ymax) then
    write(*,*) 'Error in map_interpolation: ymax=',grid_ymax,' and y=',y
    return
  endif

  k0=nx/2
  l0=ny/2
  status=-1
  
  do attempt=1,10

    kprev=k0
    lprev=l0

    ux = xgr(k0+1,l0)-xgr(k0,l0)
    uy = ygr(k0+1,l0)-ygr(k0,l0)

    vx = xgr(k0,l0+1)-xgr(k0,l0)
    vy = ygr(k0,l0+1)-ygr(k0,l0)

    dx=x-xgr(k0,l0)
    dy=y-ygr(k0,l0)

    d=ux*vy-uy*vx

    a=(dx*vy-dy*vx)/d
    b=(ux*dy-uy*dx)/d
    
    k0 = k0 + nint (a)
    l0 = l0 + nint (b)

    if(l0<1) l0=1
    if(l0>ny) l0=ny-1

    if(l0<1) l0=1
    if(l0>ny) l0=ny-1
    
!    write(*,*) attempt,k0,l0
    
    if((kprev==k0) .and. (lprev==l0)) exit

   enddo
   
! initial guess on the closest point

    ux = xgr(k0+1,l0)-xgr(k0,l0)
    uy = ygr(k0+1,l0)-ygr(k0,l0)

    vx = xgr(k0,l0+1)-xgr(k0,l0)
    vy = ygr(k0,l0+1)-ygr(k0,l0)

    dx=x-xgr(k0,l0)
    dy=y-ygr(k0,l0)

    d=ux*vy-uy*vx

    a=(dx*vy-dy*vx)/d
    b=(ux*dy-uy*dx)/d
    
!    write(*,*) a,b

    k0 = k0 + nint (a)
    l0 = l0 + nint (b)
    
!    write(*,*) k0,l0
  
  
! recompute the a,b in first sector

    ux = xgr(k0+1,l0)-xgr(k0,l0)
    uy = ygr(k0+1,l0)-ygr(k0,l0)

    vx = xgr(k0,l0+1)-xgr(k0,l0)
    vy = ygr(k0,l0+1)-ygr(k0,l0)

    dx=x-xgr(k0,l0)
    dy=y-ygr(k0,l0)

    d=ux*vy-uy*vx

    a=(dx*vy-dy*vx)/d
    b=(ux*dy-uy*dx)/d
  
!  write(*,*) '1 sec',a,b

  step_k=1
  step_l=1

  if (a>-eps .and. b>-eps) goto 1000

! recompute the a,b in second sector

    ux = xgr(k0-1,l0)-xgr(k0,l0)
    uy = ygr(k0-1,l0)-ygr(k0,l0)

    vx = xgr(k0,l0+1)-xgr(k0,l0)
    vy = ygr(k0,l0+1)-ygr(k0,l0)

    dx=x-xgr(k0,l0)
    dy=y-ygr(k0,l0)

    d=ux*vy-uy*vx

    a=(dx*vy-dy*vx)/d
    b=(ux*dy-uy*dx)/d
  
!  write(*,*) '2 sec',a,b

  step_k=-1
  step_l=1

  if (a>-eps .and. b>-eps) goto 1000

! recompute the a,b in third sector

    ux = xgr(k0-1,l0)-xgr(k0,l0)
    uy = ygr(k0-1,l0)-ygr(k0,l0)

    vx = xgr(k0,l0-1)-xgr(k0,l0)
    vy = ygr(k0,l0-1)-ygr(k0,l0)

    dx=x-xgr(k0,l0)
    dy=y-ygr(k0,l0)

    d=ux*vy-uy*vx

    a=(dx*vy-dy*vx)/d
    b=(ux*dy-uy*dx)/d
  
!  write(*,*) '3 sec',a,b

  step_k=-1
  step_l=-1

  if (a>-eps .and. b>-eps) goto 1000

! recompute the a,b in fourth sector

    ux = xgr(k0+1,l0)-xgr(k0,l0)
    uy = ygr(k0+1,l0)-ygr(k0,l0)

    vx = xgr(k0,l0-1)-xgr(k0,l0)
    vy = ygr(k0,l0-1)-ygr(k0,l0)

    dx=x-xgr(k0,l0)
    dy=y-ygr(k0,l0)

    d=ux*vy-uy*vx

    a=(dx*vy-dy*vx)/d
    b=(ux*dy-uy*dx)/d
  
!  write(*,*) '4 sec',a,b

  step_k=1
  step_l=-1

  if (a>-eps .and. b>-eps) then
   goto 1000
  else
   write(*,*) 'issues with interpolation in sectors'
   return
  endif

!-----------------------------------------------------------------------
! Interpolation (bilinear)

1000 continue

  count=0

  xx(0)=1.-a
  xx(1)=a
    
  yy(0)=1.-b
  yy(1)=b

  count=0
  do i=0,1
    do j=0,1
      count = count + 1
      i1=step_k*i+k0; j1=step_l*j+l0
      indexij(count,i2,j2)=i1+nx*(j1-1)
      weight(count,i2,j2)=xx(i)*yy(j)
    enddo
  enddo

  w=0.d0
  do m=1,count
    w=w+weight(m,i2,j2)
  enddo

  if (w > 0.d0) then
    do m=1,count
      weight(m,i2,j2)=weight(m,i2,j2)/w
    enddo
  else
    do m=1,count
      weight(m,i2,j2)=mask
    enddo
  endif

 enddo
enddo

 return
end subroutine weight_interpolation_irreg


END MODULE map
