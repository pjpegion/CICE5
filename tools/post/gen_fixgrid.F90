program gen_fixgrid
!
! Denise.Worthen@noaa.gov
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this is a stripped down and modified version of the code used 
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
! Area of the T-grid cell is obtained as in MOM_grid_initialize where
! tmpV = dx on SG and tmpU is dy on SG
!
!    dxT(i,j) = tmpV(i2-1,j2-1) + tmpV(i2,j2-1)
!    dyT(i,j) = tmpU(i2-1,j2-1) + tmpU(i2-1,j2)
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
! Vertices are defined counter-clockwise from upper right. Ct-grid vertices
! are located on the Bu grid; Cu vertices on the Cv grid, Cv vertices on the Cu
! grid and Bu vertices on the Ct grid. For example, for the Ct-grid, the vertices
! are: 
!             Vertex #2             Vertex #1
!             Bu(i-1,j)             Bu(i,j)  
!                         Ct(i,j)
!           Bu(i-1,j-1)             Bu(i,j-1)
!             Vertex #3             Vertex #4
!
! so that the vertices of any Ct(i,j) are found as off-sets of the i,j index on the
! Bu grid 
! 
!     iVertCt(4) = (/0, -1, -1, 0/)
!     jVertCt(4) = (/0, 0, -1, -1/)
! 
! Careful examination of the Cu,Cv and Bu grids lead to similar definitions for the
! i,j offsets required to extract the other grid stragger vertices locations, all of
! which can be defined in terms of the iVertCt and jVertCt values
!  
! Special treatment is require at the bottom of the grid, where the verticies of the
! Ctand Cu grid must be set manually (note, these points are on land.) The top of 
! the grid also requires special treatment because the required verticies are located
! across the tripole seam. This is accomplished by creating 1-d arrays which hold
! the Ct and Cu grid point locations on the matched seam.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use param
  use grdvar
  use debugprint
  use fixgriddefs
  use netcdf

  implicit none

  character(len=256) :: dirsrc = '/scratch4/NCEPDEV/nems/noscrub/emc.nemspara/RT/FV3-MOM6-CICE5/benchmark-20180913/MOM6_FIX_025deg/'
  character(len=256) :: dirout = '/scratch4/NCEPDEV/ocean/save/Denise.Worthen/NEMS_INPUT0.1/ocnicepost/'
  character(len= 10) :: res = 'mx025'

  ! super-grid source variables
  real(kind=8), dimension(0:nx,0:ny)   ::  x,  y
  real(kind=8), dimension(  nx,0:ny)   :: dx
  real(kind=8), dimension(0:nx,  ny)   :: dy

  real(kind=8) :: dxT, dyT

  character(len=256) :: fname_out, fname_in
  character(len=256) :: history
  character(len=  8) :: cdate

  integer :: rc,ncid
  integer :: id, dim2(2), dim3(3)
  integer :: ni_dim,nj_dim,nv_dim
  integer :: i,j,n,ii,jj,i2,j2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up the arrays to retrieve the vertices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  iVertCu = iVertCt + 1; jVertCu = jVertCt + 0
  iVertCv = iVertCt + 0; jVertCv = jVertCt + 1
  iVertBu = iVertCt + 1; jVertBu = jVertCt + 1

  print '(a8,4i6)','iVertCt ',(iVertCt(i),i=1,4)
  print '(a8,4i6)','jVertCt ',(jVertCt(i),i=1,4)
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

  latCt_vert = -9999.0 ; lonCt_vert = -9999.0
  latCu_vert = -9999.0 ; lonCu_vert = -9999.0
  latCv_vert = -9999.0 ; lonCv_vert = -9999.0
  latBu_vert = -9999.0 ; lonBu_vert = -9999.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the land mask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  fname_in = trim(dirsrc)//"ocean_mask.nc"
  rc = nf90_open(fname_in, nf90_nowrite, ncid)
  rc = nf90_inq_varid(ncid, 'mask',    id) 
  rc = nf90_get_var(ncid,       id, latCt) !temp use
  rc = nf90_close(ncid)

  wet = int(latCt,4)

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

  rc = nf90_inq_varid(ncid, 'dx', id) 
  rc = nf90_get_var(ncid,     id, dx)

  rc = nf90_inq_varid(ncid, 'dy', id) 
  rc = nf90_get_var(ncid,     id, dy)

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
     lonCt(i,j) =     x(i2-1,j2-1)
     lonCu(i,j) =     x(i2,  j2-1)
     lonCv(i,j) =     x(i2-1,j2  )
    !deg
     latCt(i,j) =     y(i2-1,j2-1)
     latCu(i,j) =     y(i2,  j2-1)
     latCv(i,j) =     y(i2-1,j2  )
    !m2
            dxT = dx(i2-1,j2-1) + dx(i2,j2-1)
            dyT = dy(i2-1,j2-1) + dy(i2-1,j2)
    areaCt(i,j) = dxT*dyT
   enddo
  enddo

  where(lonCt .lt. 0.0)lonCt = lonCt + 360.d0
  where(lonCu .lt. 0.0)lonCu = lonCu + 360.d0
  where(lonCv .lt. 0.0)lonCv = lonCv + 360.d0
  where(lonBu .lt. 0.0)lonBu = lonBu + 360.d0
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

  call checkseam

  do i = 1,ni
    i2 = ipole(2)+(ipole(1)-i)+1
    xlonCt(i) = lonCt(i2,nj)
    xlatCt(i) = latCt(i2,nj)
  enddo

  !do i = 1,10
  !  i2 = ipole(2)+(ipole(1)-i)
  !  print *,i,i2,lonCu(i,nj)
  !enddo
  !do i = 1430,1440
  !  i2 = ipole(2)+(ipole(1)-i)
  !  if(i2 .lt. 1)i2 = ni
  !  print *,i,i2,lonCu(i,nj)
  !enddo
 
  do i = 1,ni
    i2 = ipole(2)+(ipole(1)-i)
    if(i2 .lt. 1)i2 = ni
   xlonCu(i) = lonCu(i2,nj)
   xlatCu(i) = latCu(i2,nj)
   !print *,i,xlonCu(i),lonCu(i2,nj)
  enddo
 
  call checkxlatlon
  !do i = 1,ni
  !  i2 = ipole(2)+(ipole(1)-i)
  !  if(i2 .lt. 1)i2 = ni
  ! print *,i,xlonCu(i),lonCu(i2,nj)
  !enddo

  !approx lat at grid bottom 
  do i = 1,ni
   dlatBu(i) = latBu(i,1) + 2.0*(latCu(i,1) - latBu(i,1))
   dlatCv(i) = latCt(i,1) + 2.0*(latCt(i,1) - latCv(i,1))
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fill grid vertices variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Ct and Cu grids align in j 
  call fill_vertices(2,nj  , iVertCt,jVertCt, latBu,lonBu, latCt_vert,lonCt_vert)
  call           fill_bottom(iVertCt,jVertCt, latBu,lonBu, latCt_vert,lonCt_vert,dlatBu)

  call fill_vertices(2,nj  , iVertCu,jVertCu, latCv,lonCv, latCu_vert,lonCu_vert)
  call           fill_bottom(iVertCu,jVertCu, latCv,lonCv, latCu_vert,lonCu_vert,dlatCv)

  !Cv and Bu grids align in j
  call fill_vertices(1,nj-1, iVertCv,jVertCv, latCu,lonCu, latCv_vert,lonCv_vert)
  call              fill_top(iVertCv,jVertCv, latCu,lonCu, latCv_vert,lonCv_vert, xlatCu, xlonCu)

  call fill_vertices(1,nj-1, iVertBu,jVertBu, latCt,lonCt, latBu_vert,lonBu_vert)
  call              fill_top(iVertBu,jVertBu, latCt,lonCt, latBu_vert,lonBu_vert, xlatCt, xlonCt)
 
  call checkpoint

  if(minval(latCt_vert) .lt. -1.e3)stop
  if(minval(lonCt_vert) .lt. -1.e3)stop
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
  fname_out= trim(dirout)//'tripole.'//trim(res)//'.nc'
  print *,trim(fname_out)

  ! create the file
  rc = nf90_create(trim(fname_out), nf90_write, ncid)
  print *, 'writing grid to ',trim(fname_out)
  print *, 'nf90_create = ',trim(nf90_strerror(rc))

  rc = nf90_def_dim(ncid,'ni', ni, ni_dim)
  rc = nf90_def_dim(ncid,'nj', nj, nj_dim)
  rc = nf90_def_dim(ncid,'nv', nv, nv_dim)
  
  !mask
  dim2(2) = nj_dim
  dim2(1) = ni_dim
   rc = nf90_def_var(ncid, 'wet', nf90_int, dim2, id)

  !area
  dim2(2) = nj_dim
  dim2(1) = ni_dim
   rc = nf90_def_var(ncid, 'area', nf90_double, dim2, id)
   rc = nf90_put_att(ncid, id,     'units',  'm2')

  dim2(2) = nj_dim
  dim2(1) = ni_dim
  do ii = 1,ncoord
   rc = nf90_def_var(ncid, trim(fixgrid(ii)%var_name), nf90_double, dim2, id)
   rc = nf90_put_att(ncid, id,     'units', trim(fixgrid(ii)%unit_name))
   rc = nf90_put_att(ncid, id, 'long_name', trim(fixgrid(ii)%long_name))
   if(trim(fixgrid(ii)%var_name(1:3)) .eq. "lon")then
    rc = nf90_put_att(ncid, id,  'lon_bounds', trim(fixgrid(ii)%vertices))
   else
    rc = nf90_put_att(ncid, id,  'lat_bounds', trim(fixgrid(ii)%vertices))
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

  rc = nf90_inq_varid(ncid,   'wet',      id)
  rc = nf90_put_var(ncid,        id,     wet)

  rc = nf90_inq_varid(ncid,  'area',      id)
  rc = nf90_put_var(ncid,        id,  areaCt)

  rc = nf90_inq_varid(ncid,  'lonCt',     id)
  rc = nf90_put_var(ncid,        id,   lonCt)

  rc = nf90_inq_varid(ncid,  'latCt',     id)
  rc = nf90_put_var(ncid,        id,   latCt)

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
  rc = nf90_put_var(ncid,         id,  lonCt_vert)

  rc = nf90_inq_varid(ncid,  'latCt_vert',     id)
  rc = nf90_put_var(ncid,         id,  latCt_vert)

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
