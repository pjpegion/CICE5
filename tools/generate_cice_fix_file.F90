program generate_cice_fix_file

! Denise.Worthen@noaa.gov
! Jiande.Wang@noaa.gov
!
!#define output_grid_1deg
#define output_grid_qdeg
! writes out additional variables not needed by CICE but which can be 
! helpful in diagnosing grid generation errors
#define debug
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this code generate CICE gird fixed file based on MOM6 ocean_hgrid.nc
! information on MOM6 supergrid can be found at
! https://gist.github.com/adcroft/c1e207024fe1189b43dddc5f1fe7dd6c
!
! also: https://mom6.readthedocs.io/en/latest/api/generated/modules/mom_grid.html
!
! also:
! MOM_grid_initialize.F90 :
!  MOM6 variable geoLonBu <==> CICE variable ulon
!  MOM6 variable geoLatBu <==> CICE variable ulat
!  MOM6 variable     dxCv <==> CICE variable htn
!  MOM6 variable     dyCu <==> CICE variable hte
!
! MOM6 code snippets follow:
!
! from MOM_grid_initialize.F90  (tmpZ = x)
!  do J=G%JsdB,G%JedB ; do I=G%IsdB,G%IedB ; i2 = 2*I ; j2 = 2*J
!    G%geoLonBu(I,J) = tmpZ(i2,j2)
! so....
!          ulon(I,J) = x(i2,j2)
 
! from MOM_grid_initialize.F90  (tmpZ = y)
!  do J=G%JsdB,G%JedB ; do I=G%IsdB,G%IedB ; i2 = 2*I ; j2 = 2*J
!    G%geoLatBu(I,J) = tmpZ(i2,j2)
! so....
!          ulat(I,J) = y(i2,j2)

! from MOM_grid_initialize.F90  (tmpV = dx)
!  do J=G%JsdB,G%JedB ; do i=G%isd,G%ied ; i2 = 2*i ; j2 = 2*j
!    dxCv(i,J) = tmpV(i2-1,j2) + tmpV(i2,j2)
! so....
!     htn(i,J) =   dx(i2-1,j2) +   dx(i2,j2)

! from MOM_grid_initialize.F90  (tmpU = dy)
!  do J=G%JsdB,G%JedB ; do i=G%isd,G%ied ; i2 = 2*i ; j2 = 2*j
!    dyCu(I,j) = tmpU(i2,j2-1) + tmpU(i2,j2)
! so....
!     hte(I,j) =   dy(i2,j2-1) +   dy(i2,j2)
!
! rotation angle on supergrid vertices can be found 
! using the formula in MOM_shared_initialization.F90, accounting 
! for indexing difference between reduced grid and super grid
!
!         SuperGrid                 Reduced grid
! 
!  i-1,j+1         i+1,j+1
!     X-------X-------X             I-1,J      I,J
!     |       |       |                X-------X      
!     |       |       |                |       |
!     |       | i,j   |                |   T   |
!     X-------X-------X                |       |
!     |       |       |                X-------X
!     |       |       |             I-1,J-1   I,J-1
!     |       |       |
!     X-------X-------X
!  i-1,j-1         i+1,j-1
!      
! so that in angle formulae 
!         I==>i+1,I-1==>i-1
!         J==>j+1,J-1==>j-1 
!
! CICE expects angle to be XY -> LatLon so change the sign from MOM6 value
! This has been determined from the HYCOM code: ALL/cice/src/grid2cice.f 
!
!            anglet(i,j) =    -pang(i+i0,  j+j0)   !radians
!c           pang is from lon-lat to x-y, but anglet is the reverse
!
! where anglet is the angle variable being written to the CICE grid file 
! and pang is HYCOM's own rotation angle. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use netcdf
  use icegriddefs

  implicit none
#ifdef output_grid_1deg
  integer, parameter :: ni = 360, nj = 210
  character(len=256) :: dirsrc = '/scratch3/NCEPDEV/stmp2/Jiande.Wang/UFS-1x1/master-20180521_1x1/MOM6_FIX_1deg/'
  character(len= 10) :: dirout = './'
  character(len= 10) :: res = 'mx1'
#endif
#ifdef output_grid_qdeg
  integer, parameter :: ni = 1440, nj = 1080
  character(len=256) :: dirsrc = '/scratch4/NCEPDEV/nems/noscrub/emc.nemspara/RT/FV3-MOM6-CICE5/master-20180821/MOM6_FIX_025deg/'
  character(len=256) :: dirout = './'
  character(len= 10) :: res = 'mx025'
#endif

  real(kind=8), parameter ::      pi = 4.0*atan(1.0)
  real(kind=8), parameter :: deg2rad = pi/180.0d0
  
  ! super-grid source variables
  integer, parameter :: nx  = ni*2, ny  = nj*2
  real(kind=8), dimension(0:nx,0:ny)   :: x, y, angq 
  real(kind=8), dimension(  nx,0:ny)   :: dx
  real(kind=8), dimension(0:nx,  ny)   :: dy

  ! required CICE grid variables
  real(kind=8), dimension(ni,nj) :: ulon, ulat 
  real(kind=8), dimension(ni,nj) :: dxt, dyt, htn, hte
  real(kind=8), dimension(ni,nj) :: angle
#ifdef debug
  real(kind=8), dimension(ni,nj) ::  latT, lonT  ! lat and lon of T on C-grid
  real(kind=8), dimension(ni,nj) :: latCv, lonCv ! lat and lon of V on C-grid
  real(kind=8), dimension(ni,nj) :: latCu, lonCu ! lat and lon of U on C-grid
  real(kind=8), dimension(ni,nj) :: anglet
#endif

  integer, parameter :: ncice = 5  & ! required
#ifdef debug
                              + 7    ! extra
#else
                              + 0    ! extra
#endif
  character(len=256) :: fname_out, fname_in

  real(kind=8) :: lon_scale
  integer :: status,ncid
  integer :: id, vardim(2)
  integer :: ni_dim,nj_dim
  integer :: i,j,ii,jj,i2,j2
  integer :: ipole(2),i1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read supergrid file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  fname_in = trim(dirsrc)//'ocean_hgrid.nc'

  status = nf90_open(fname_in, nf90_nowrite, ncid)
  print *, 'reading supergrid from ',trim(fname_in)
  print *, 'nf90_open = ',trim(nf90_strerror(status))
  
  status = nf90_inq_varid(ncid, 'x', id)  !lon
  status = nf90_get_var(ncid,    id,  x)
 
  status = nf90_inq_varid(ncid, 'y', id)  !lat
  status = nf90_get_var(ncid,    id,  y)
 
  status = nf90_inq_varid(ncid, 'dx', id)  
  status = nf90_get_var(ncid,     id, dx)
  
  status = nf90_inq_varid(ncid, 'dy', id)  
  status = nf90_get_var(ncid,     id, dy)
  status = nf90_close(ncid)

  do j=1,ny-1 ; do i=1,nx-1
    lon_scale    = cos((y(i-1,j-1) + y(i+1,j-1  ) + &
                        y(i-1,j+1) + y(i+1,j+1)) * atan(1.0)/180)
    angq(i,j)    = atan2((x(i-1,j+1) + x(i+1,j+1) - &
                          x(i-1,j-1) - x(i+1,j-1))*lon_scale, &
                          y(i-1,j+1) + y(i+1,j+1) - &
                          y(i-1,j-1) - y(i+1,j-1) )
  enddo ; enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fill cice grid variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do j = 1,nj
   do i = 1,ni
     i2 = 2*i ; j2 = 2*j
    !deg->rad
     ulon(i,j) =     x(i2,j2)*deg2rad
     ulat(i,j) =     y(i2,j2)*deg2rad
    !in rad already
    angle(i,j) = -angq(i2,j2)
    !m->cm
      htn(i,j) = (dx(i2-1,j2) + dx(i2,j2))*100.0
      hte(i,j) = (dy(i2,j2-1) + dy(i2,j2))*100.0
#ifdef debug
    !deg
      lonT(i,j) =     x(i2-1,j2-1)
     lonCu(i,j) =     x(i2,  j2-1)
     lonCv(i,j) =     x(i2-1,j2  )
    !deg
      latT(i,j) =     y(i2-1,j2-1)
     latCu(i,j) =     y(i2,  j2-1)
     latCv(i,j) =     y(i2-1,j2  )
    !in rad already
    anglet(i,j) = -angq(i2-1,j2-1)
#endif
   enddo
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For the 1/4deg grid, hte at j=720 and j = 1440 is identically=0.0 for 
! j > 840 (64.0N). These are land points, but since CICE uses hte to 
! generate remaining variables, setting them to zero will cause problems
! For 1deg grid, hte at ni/2 and ni are very small O~10-12, so test for 
! hte < 1.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  print *,'min vals of hte at folds ',minval(hte(ni/2,:)),minval(hte(ni,:))
  do j = 1,nj
     ii = ni/2
   if(hte(ii,j) .le. 1.0)hte(ii,j) = 0.5*(hte(ii-1,j) + hte(ii+1,j))
     ii = ni
   if(hte(ii,j) .le. 1.0)hte(ii,j) = 0.5*(hte(ii-1,j) + hte(   1,j))
  enddo
  print *,'min vals of hte at folds ',minval(hte(ni/2,:)),minval(hte(ni,:))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! some basic error checking
! find the i-th index of the poles at j= nj
! the corner points must lie on the pole
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ipole = -1
      j = nj
  do i = 1,ni/2
   if(ulat(i,j)/deg2rad .eq. 90.0)ipole(1) = i
  enddo
  do i = ni/2+1,ni
   if(ulat(i,j)/deg2rad .eq. 90.0)ipole(2) = i
  enddo
  print *,'poles found at ',ipole

  !htn must be the same along seam
   j = nj
  i1 = ipole(1); i2 = ipole(2)+1
  print *,'HTN across seam '
  print *,htn(i1-2,j),htn(i2+2,j)
  print *,htn(i1-1,j),htn(i2+1,j)
  print *,htn(i1,  j),htn(i2,  j)
  print *,htn(i1+1,j),htn(i2-1,j)
  print *,htn(i1+2,j),htn(i2-2,j)
#ifdef debug
  print *,'latCv across seam '
  print *,latCv(i1-2,j),latCv(i2+2,j)
  print *,latCv(i1-1,j),latCv(i2+1,j)
  print *,latCv(i1,  j),latCv(i2,  j)
  print *,latCv(i1+1,j),latCv(i2-1,j)
  print *,latCv(i1+2,j),latCv(i2-2,j)

  print *,'lonCv across seam '
  print *,lonCv(i1-2,j),lonCv(i2+2,j)
  print *,lonCv(i1-1,j),lonCv(i2+1,j)
  print *,lonCv(i1,  j),lonCv(i2,  j)
  print *,lonCv(i1+1,j),lonCv(i2-1,j)
  print *,lonCv(i1+2,j),lonCv(i2-2,j)

  print *,'anglet across seam '
  print *,angleT(i1-2,j),angleT(i2+2,j)
  print *,angleT(i1-1,j),angleT(i2+1,j)
  print *,angleT(i1,  j),angleT(i2,  j)
  print *,angleT(i1+1,j),angleT(i2-1,j)
  print *,angleT(i1+2,j),angleT(i2-2,j)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write out cice grid file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call ice_typedefine

  fname_out= trim(dirout)//'grid_cice_NEMS_'//trim(res)//'.nc'

  status = nf90_create(fname_out, nf90_write, ncid)
  print *, 'writing CICE grid to ',trim(fname_out)
  print *, 'nf90_create = ',trim(nf90_strerror(status))

  status = nf90_def_dim(ncid,'ni', ni, ni_dim)
  status = nf90_def_dim(ncid,'nj', nj, nj_dim)

  vardim(2) = nj_dim
  vardim(1) = ni_dim
  do ii = 1,ncice
   status = nf90_def_var(ncid, trim(icegrid(ii)%var_name), nf90_double, vardim, id)
   status = nf90_put_att(ncid, id,     'units', trim(icegrid(ii)%unit_name))
   status = nf90_put_att(ncid, id, 'long_name', trim(icegrid(ii)%long_name))
  enddo
   status = nf90_enddef(ncid)

  status = nf90_inq_varid(ncid,  'ulon',      id)
  status = nf90_put_var(ncid,        id,    ulon)

  status = nf90_inq_varid(ncid,  'ulat',      id)
  status = nf90_put_var(ncid,        id,    ulat)

  status = nf90_inq_varid(ncid,   'htn',      id)
  status = nf90_put_var(ncid,        id,     htn)

  status = nf90_inq_varid(ncid,   'hte',      id)
  status = nf90_put_var(ncid,        id,     hte)
 
  status = nf90_inq_varid(ncid,  'angle',     id)
  status = nf90_put_var(ncid,         id,  angle)
#ifdef debug
  status = nf90_inq_varid(ncid, 'anglet',     id)
  status = nf90_put_var(ncid,         id, anglet)

  status = nf90_inq_varid(ncid,  'lonT',     id)
  status = nf90_put_var(ncid,        id,   lonT)

  status = nf90_inq_varid(ncid,  'latT',     id)
  status = nf90_put_var(ncid,        id,   latT)

  status = nf90_inq_varid(ncid, 'lonCv',      id)
  status = nf90_put_var(ncid,        id,   lonCv)

  status = nf90_inq_varid(ncid, 'latCv',      id)
  status = nf90_put_var(ncid,        id,   latCv)
  
  status = nf90_inq_varid(ncid, 'lonCu',      id)
  status = nf90_put_var(ncid,        id,   lonCu)

  status = nf90_inq_varid(ncid, 'latCu',      id)
  status = nf90_put_var(ncid,        id,   latCu)
#endif
  status = nf90_close(ncid)
end program generate_cice_fix_file
