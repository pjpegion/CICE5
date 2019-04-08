program gen_fixgrid
!
! Denise.Worthen@noaa.gov
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this is a stripped down and slightly modified version of the code used 
! to generate the CICE5 grid from the MOM6 supergrid (ocean_hgrid.nc)
!
! this code generates a fixed grid file and is used to create the interpolation
! weights for ice/ocean post
!
! information on MOM6 supergrid can be found at
! https://gist.github.com/adcroft/c1e207024fe1189b43dddc5f1fe7dd6c
!
! also: https://mom6.readthedocs.io/en/latest/api/generated/modules/mom_grid.html
!
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
!
! Tripole Seam flip: ipL,ipR left,right poles on seam
!
! ipL-1     ipL    ipL+1       ipR-1     ipR    ipR+1
!    x-------x-------x     |||    x-------x-------x 
!
! Fold over; ipL must align with ipR
!  
!  ipR+1     ipR    ipR-1
!     x-------x-------x 
!  ipL-1     ipL    ipL+1
!     x-------x-------x 
!
!
! Vertices are defined:
!   
!             vert(2,i,j)             vert(1,i,j)  
!                            T(i,j)
!             vert(3,i,j)             vert(4,i,j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use netcdf
  use fixgriddefs

  implicit none
  integer, parameter :: ni = 1440, nj = 1080, nv = 4
  character(len=256) :: dirsrc = '/scratch4/NCEPDEV/nems/noscrub/emc.nemspara/RT/FV3-MOM6-CICE5/benchmark-20180913/MOM6_FIX_025deg/'
  character(len=256) :: dirout = './'
  character(len= 10) :: res = 'mx025'

  real(kind=8), parameter ::      pi = 4.0*atan(1.0)
  real(kind=8), parameter :: deg2rad = pi/180.0d0
  
  ! super-grid source variables
  integer, parameter :: nx  = ni*2, ny  = nj*2
  real(kind=8), dimension(0:nx,0:ny)   :: x, y

  real(kind=8), dimension(ni,nj) ::  latT, lonT  ! lat and lon of T on C-grid
  real(kind=8), dimension(ni,nj) :: latCv, lonCv ! lat and lon of V on C-grid
  real(kind=8), dimension(ni,nj) :: latCu, lonCu ! lat and lon of U on C-grid
  real(kind=8), dimension(ni,nj) :: latBu, lonBu ! lat and lon of corners on C-grid

  real(kind=8), dimension(ni,nj,nv) ::  latT_vert,  lonT_vert
  real(kind=8), dimension(ni,nj,nv) :: latCu_vert, lonCu_vert
  real(kind=8), dimension(ni,nj,nv) :: latCv_vert, lonCv_vert
  real(kind=8), dimension(ni,nj,nv) :: latBu_vert, lonBu_vert

  ! ij offsets moving counter-clockwise around each T(i,j)
  integer, dimension(nv) :: iVertT = (/0, -1, -1,  0/)
  integer, dimension(nv) :: jVertT = (/0,  0, -1, -1/)

  integer, dimension(nv) :: iVertBu, iVertCu, iVertCv
  integer, dimension(nv) :: jVertBu, jVertCu, jVertCv 

  ! need across seam values of T,Cu points to retrieve vertices of Bu and Cv grids
  real(kind=8), dimension(ni) ::  xlonT,  xlatT
  real(kind=8), dimension(ni) :: xlonCu, xlatCu

  integer, parameter :: ncoord = 2*4             ! 4sets of lat/lon pairs
  integer, parameter :: nverts = 2*4             ! 4sets of lat/lon pairs vertices
  integer, parameter ::  nvars = ncoord + nverts

  character(len=256) :: fname_out, fname_in
  character(len=256) :: history
  character(len=  8) :: cdate

  integer :: rc,ncid
  integer :: id, dim2(2), dim3(3)
  integer :: ni_dim,nj_dim,nv_dim
  integer :: i,j,n,m,ii,jj,i2,j2,ip1,im1,jp1,jm1
  integer :: ipole(2),i1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up the arrays to retrieve the vertices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  iVertCu = iVertT + 1; jVertCu = jVertT + 0
  iVertCv = iVertT + 0; jVertCv = jVertT + 1
  iVertBu = iVertT + 1; jVertBu = jVertT + 1

  print '(a8,4i6)','iVertT  ',( iVertT(i),i=1,4)
  print '(a8,4i6)','jVertT  ',( jVertT(i),i=1,4)
  print *
  print '(a8,4i6)','iVertCu ',(iVertCu(i),i=1,4)
  print '(a8,4i6)','jVertCu ',(jVertCu(i),i=1,4)
  print *
  print '(a8,4i6)','iVertCv ',(iVertCv(i),i=1,4)
  print '(a8,4i6)','jVertCv ',(jVertCv(i),i=1,4)
  print *
  print '(a8,4i6)','iVertBu ',(iVertBu(i),i=1,4)
  print '(a8,4i6)','jVertBu ',(jVertBu(i),i=1,4)
  print *

   latT_vert = -9999.0 ;  lonT_vert = -9999.0
  latCu_vert = -9999.0 ; lonCu_vert = -9999.0
  latCv_vert = -9999.0 ; lonCv_vert = -9999.0
  latBu_vert = -9999.0 ; lonBu_vert = -9999.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read supergrid file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  fname_in = trim(dirsrc)//'ocean_hgrid.nc'

  rc = nf90_open(fname_in, nf90_nowrite, ncid)
  print *, 'reading supergrid from ',trim(fname_in)
  print *, 'nf90_open = ',trim(nf90_strerror(rc))
  
  rc = nf90_inq_varid(ncid, 'x', id)  !lon
  rc = nf90_get_var(ncid,    id,  x)
 
  rc = nf90_inq_varid(ncid, 'y', id)  !lat
  rc = nf90_get_var(ncid,    id,  y)
  rc = nf90_close(ncid)
  print *,'super grid size ',size(y,1),size(y,2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fill grid variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do j = 1,nj
   do i = 1,ni
     i2 = 2*i ; j2 = 2*j
    !deg
     lonBu(i,j) =     x(i2,j2)
     latBu(i,j) =     y(i2,j2)
    !deg
      lonT(i,j) =     x(i2-1,j2-1)
     lonCu(i,j) =     x(i2,  j2-1)
     lonCv(i,j) =     x(i2-1,j2  )
    !deg
      latT(i,j) =     y(i2-1,j2-1)
     latCu(i,j) =     y(i2,  j2-1)
     latCv(i,j) =     y(i2-1,j2  )
   enddo
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! some basic error checking
! find the i-th index of the poles at j= nj
! the corner points must lie on the pole
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ipole = -1
      j = nj
  do i = 1,ni/2
   if(latBu(i,j) .eq. 90.0)ipole(1) = i
  enddo
  do i = ni/2+1,ni
   if(latBu(i,j) .eq. 90.0)ipole(2) = i
  enddo
  print *,'poles found at ',ipole,latBu(ipole(1),nj),latBu(ipole(2),nj)

  j = nj
  i1 = ipole(1); i2 = ipole(2)+1
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

  print *,'latCu across seam '
  print *,latCu(i1-2,j),latCu(i2+1,j)
  print *,latCu(i1-1,j),latCu(i2+0,j)
  print *,latCu(i1,  j),latCu(i2-1,j)
  print *,latCu(i1+1,j),latCu(i2-2,j)
  print *,latCu(i1+2,j),latCu(i2-3,j)

  print *,'lonCu across seam '
  print *,lonCu(i1-2,j),lonCu(i2+1,j)
  print *,lonCu(i1-1,j),lonCu(i2+0,j)
  print *,lonCu(i1,  j),lonCu(i2-1,j)
  print *,lonCu(i1+1,j),lonCu(i2-2,j)
  print *,lonCu(i1+2,j),lonCu(i2-3,j)

  print *,'latT across seam '
  print *,latT(i1-2,j),latT(i2+2,j)
  print *,latT(i1-1,j),latT(i2+1,j)
  print *,latT(i1,  j),latT(i2,  j)
  print *,latT(i1+1,j),latT(i2-1,j)
  print *,latT(i1+2,j),latT(i2-2,j)

  print *,'lonT across seam '
  print *,lonT(i1-2,j),lonT(i2+2,j)
  print *,lonT(i1-1,j),lonT(i2+1,j)
  print *,lonT(i1,  j),lonT(i2,  j)
  print *,lonT(i1+1,j),lonT(i2-1,j)
  print *,lonT(i1+2,j),lonT(i2-2,j)
  print *

  do i = 1,ni
    i2 = ipole(2)+(ipole(1)-i)+1
    xlonT(i) = lonT(i2,nj)
    xlatT(i) = latT(i2,nj)
  enddo

  do i = 1,ni
    i2 = ipole(2)+(ipole(1)-i)
    if(i2 .lt. 1)i2 = ni
   xlonCu(i) = lonCu(i2,nj)
   xlatCu(i) = latCu(i2,nj)
  enddo
  
  print *,'============== T grid ==============='
  do i = 355,365
   print '(i5,6f12.5)',i,lonT(i,nj),xlonT(i),lonT(i,nj)+xlonT(i),latT(i,nj),xlatT(i),latT(i,nj)-xlatT(i)
  enddo
  print *
  do i = 355,365
   print '(i5,6f12.5)',i,lonT(i,nj),xlonT(i),lonT(i,nj)+xlonT(i),latT(i,nj),xlatT(i),latT(i,nj)-xlatT(i)
  enddo
  print *

  do i = 1075,1085
   print '(i5,6f12.5)',i,lonT(i,nj),xlonT(i),lonT(i,nj)+xlonT(i),latT(i,nj),xlatT(i),latT(i,nj)-xlatT(i)
  enddo
  print *
  do i = 1075,1085
   print '(i5,6f12.5)',i,lonT(i,nj),xlonT(i),lonT(i,nj)+xlonT(i),latT(i,nj),xlatT(i),latT(i,nj)-xlatT(i)
  enddo
  print *

  print *,'============== T grid ==============='
  do i = 1430,1440
   print '(i5,6f12.5)',i,lonT(i,nj),xlonT(i),lonT(i,nj)+xlonT(i),latT(i,nj),xlatT(i),latT(i,nj)-xlatT(i)
  enddo
  print *
  do i = 1430,1440
   print '(i5,6f12.5)',i,lonT(i,nj),xlonT(i),lonT(i,nj)+xlonT(i),latT(i,nj),xlatT(i),latT(i,nj)-xlatT(i)
  enddo

  print *,'============== Cu grid ==============='
  do i = 355,365
   print '(i5,6f12.5)',i,lonCu(i,nj),xlonCu(i),lonCu(i,nj)+xlonCu(i),latCu(i,nj),xlatCu(i),latCu(i,nj)-xlatCu(i)
  enddo
  print *
  do i = 355,365
   print '(i5,6f12.5)',i,lonCu(i,nj),xlonCu(i),lonCu(i,nj)+xlonCu(i),latCu(i,nj),xlatCu(i),latCu(i,nj)-xlatCu(i)
  enddo
  print *

  do i = 1075,1085
   print '(i5,6f12.5)',i,lonCu(i,nj),xlonCu(i),lonCu(i,nj)+xlonCu(i),latCu(i,nj),xlatCu(i),latCu(i,nj)-xlatCu(i)
  enddo
  print *
  do i = 1075,1085
   print '(i5,6f12.5)',i,lonCu(i,nj),xlonCu(i),lonCu(i,nj)+xlonCu(i),latCu(i,nj),xlatCu(i),latCu(i,nj)-xlatCu(i)
  enddo
  print * 

  print *,'============== Cu grid ==============='
  do i = 1430,1440
   print '(i5,6f12.5)',i,lonCu(i,nj),xlonCu(i),lonCu(i,nj)+xlonCu(i),latCu(i,nj),xlatCu(i),latCu(i,nj)-xlatCu(i)
  enddo
  print *
  do i = 1430,1440
   print '(i5,6f12.5)',i,lonCu(i,nj),xlonCu(i),lonCu(i,nj)+xlonCu(i),latCu(i,nj),xlatCu(i),latCu(i,nj)-xlatCu(i)
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fill grid vertices variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! T and Cu grid share alignment in j
  do j = 2,nj
   do i = 1,ni

    do n = 1,nv
      ii = i + iVertT(n); jj = j + jVertT(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latT_vert(i,j,n)   = latBu(ii,jj) 
      lonT_vert(i,j,n)   = lonBu(ii,jj) 
    enddo

    do n = 1,nv
      ii = i + iVertCu(n); jj = j + jVertCu(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latCu_vert(i,j,n)  = latCv(ii,jj)  
      lonCu_vert(i,j,n)  = lonCv(ii,jj)   
    enddo
   enddo
  enddo

  ! grid bottom (j=1) 
  ! vertices 1,2 are available
  ! vertices 3,4 must be set manually
  j = 1
  do i = 1,ni
    do n = 1,2
      ii = i + iVertT(n); jj = j + jVertT(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latT_vert(i,j,n)   = latBu(ii,jj)
      lonT_vert(i,j,n)   = lonBu(ii,jj)
    enddo
      lonT_vert(i,j, 3) = lonT_vert(i,j,2)
      lonT_vert(i,j, 4) = lonT_vert(i,j,1)
      latT_vert(i,j, 3) = -90.0d0
      latT_vert(i,j, 4) = -90.0d0

    do n = 1,2
      ii = i + iVertCu(n); jj = j + jVertCu(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latCu_vert(i,j,n)  = latCv(ii,jj)
      lonCu_vert(i,j,n)  = lonCv(ii,jj)
    enddo
      lonCu_vert(i,j, 3) = lonCu_vert(i,j,2)
      lonCu_vert(i,j, 4) = lonCu_vert(i,j,1)
      latCu_vert(i,j, 3) = -90.0d0
      latCu_vert(i,j, 4) = -90.0d0
  enddo !i

  ! Bu and Cv grid share alignment in j
  do j = 1,nj-1
   do i = 1,ni
    do n = 1,nv
      ii = i + iVertCv(n); jj = j + jVertCv(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latCv_vert(i,j,n)  = latCu(ii,jj)  
      lonCv_vert(i,j,n)  = lonCu(ii,jj)   
    enddo

    do n = 1,nv
      ii = i + iVertBu(n); jj = j + jVertBu(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latBu_vert(i,j,n)  =  latT(ii,jj)
      lonBu_vert(i,j,n)  =  lonT(ii,jj)
    enddo
   enddo
  enddo

  !grid top (j=nj)
  ! vertices 3,4 are available
  ! vertices 1,2 must be set manually using 'across seam' values
  j = nj
  do i = 1,ni
    do n = 3,4
      ii = i + iVertCv(n); jj = j + jVertCv(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latCv_vert(i,j,n)  = latCu(ii,jj)
      lonCv_vert(i,j,n)  = lonCu(ii,jj)
    enddo
    do n = 1,2
      ii = i + iVertCv(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latCv_vert(i,j,n)  = xlatCu(ii)
      lonCv_vert(i,j,n)  = xlonCu(ii)
    enddo

    do n = 3,4
      ii = i + iVertBu(n); jj = j + jVertBu(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latBu_vert(i,j,n)  =  latT(ii,jj)
      lonBu_vert(i,j,n)  =  lonT(ii,jj)
    enddo
    do n = 1,2
      ii = i + iVertBu(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latBu_vert(i,j,n)  = xlatT(ii)
      lonBu_vert(i,j,n)  = xlonT(ii)
    enddo
  enddo

  if(minval( latT_vert) .lt. -1.e3)stop
  if(minval( lonT_vert) .lt. -1.e3)stop
  if(minval(latCu_vert) .lt. -1.e3)stop
  if(minval(lonCu_vert) .lt. -1.e3)stop
  if(minval(latCv_vert) .lt. -1.e3)stop
  if(minval(lonCv_vert) .lt. -1.e3)stop
  if(minval(latBu_vert) .lt. -1.e3)stop
  if(minval(lonBu_vert) .lt. -1.e3)stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write out grid file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! create a history attribute
   call date_and_time(date=cdate)
   history = 'created on '//trim(cdate)//' from '//trim(fname_in)

  ! define the output variables and file name
  call fixgrid_typedefine
  fname_out= trim(dirout)//'tripole.grid.'//trim(res)//'.nc'
  print *,trim(fname_out)

  ! create the file
  rc = nf90_create(trim(fname_out), nf90_write, ncid)
  print *, 'writing grid to ',trim(fname_out)
  print *, 'nf90_create = ',trim(nf90_strerror(rc))

  rc = nf90_def_dim(ncid,'ni', ni, ni_dim)
  rc = nf90_def_dim(ncid,'nj', nj, nj_dim)
  rc = nf90_def_dim(ncid,'nv', nv, nv_dim)

  dim2(2) = nj_dim
  dim2(1) = ni_dim
  do ii = 1,ncoord
   rc = nf90_def_var(ncid, trim(fixgrid(ii)%var_name), nf90_double, dim2, id)
   rc = nf90_put_att(ncid, id,     'units', trim(fixgrid(ii)%unit_name))
   rc = nf90_put_att(ncid, id, 'long_name', trim(fixgrid(ii)%long_name))
   if(trim(fixgrid(ii)%var_name(1:3)) .eq. "lon")then
    rc = nf90_put_att(ncid, id,  'lon_bnds', trim(fixgrid(ii)%vertices))
   else
    rc = nf90_put_att(ncid, id,  'lat_bnds', trim(fixgrid(ii)%vertices))
   endif
  enddo
  dim3(3) = nv_dim
  dim3(2) = nj_dim
  dim3(1) = ni_dim
  do ii = ncoord+1,ncoord+nverts
   rc = nf90_def_var(ncid, trim(fixgrid(ii)%var_name), nf90_double, dim3, id)
   rc = nf90_put_att(ncid, id,     'units', trim(fixgrid(ii)%unit_name))
   rc = nf90_put_att(ncid, id, 'long_name', trim(fixgrid(ii)%long_name))
  enddo

   rc = nf90_put_att(ncid, nf90_global, 'history', trim(history))
   rc = nf90_enddef(ncid)

  rc = nf90_inq_varid(ncid,  'lonCt',     id)
  rc = nf90_put_var(ncid,        id,    lonT)

  rc = nf90_inq_varid(ncid,  'latCt',     id)
  rc = nf90_put_var(ncid,        id,    latT)

  rc = nf90_inq_varid(ncid, 'lonCv',      id)
  rc = nf90_put_var(ncid,        id,   lonCv)

  rc = nf90_inq_varid(ncid, 'latCv',      id)
  rc = nf90_put_var(ncid,        id,   latCv)
  
  rc = nf90_inq_varid(ncid, 'lonCu',      id)
  rc = nf90_put_var(ncid,        id,   lonCu)

  rc = nf90_inq_varid(ncid, 'latCu',      id)
  rc = nf90_put_var(ncid,        id,   latCu)

  rc = nf90_inq_varid(ncid, 'lonBu',      id)
  rc = nf90_put_var(ncid,        id,   lonBu)

  rc = nf90_inq_varid(ncid, 'latBu',      id)
  rc = nf90_put_var(ncid,        id,   latBu)

  ! vertices
  rc = nf90_inq_varid(ncid,  'lonCt_vert',     id)
  rc = nf90_put_var(ncid,         id,   lonT_vert)

  rc = nf90_inq_varid(ncid,  'latCt_vert',     id)
  rc = nf90_put_var(ncid,         id,   latT_vert)

  rc = nf90_inq_varid(ncid, 'lonCv_vert',      id)
  rc = nf90_put_var(ncid,        id,   lonCv_vert)

  rc = nf90_inq_varid(ncid, 'latCv_vert',      id)
  rc = nf90_put_var(ncid,        id,   latCv_vert)

  rc = nf90_inq_varid(ncid, 'lonCu_vert',      id)
  rc = nf90_put_var(ncid,        id,   lonCu_vert)

  rc = nf90_inq_varid(ncid, 'latCu_vert',      id)
  rc = nf90_put_var(ncid,        id,   latCu_vert)

  rc = nf90_inq_varid(ncid, 'lonBu_vert',      id)
  rc = nf90_put_var(ncid,        id,   lonBu_vert)

  rc = nf90_inq_varid(ncid, 'latBu_vert',      id)
  rc = nf90_put_var(ncid,        id,   latBu_vert)

  rc = nf90_close(ncid)

end program gen_fixgrid
