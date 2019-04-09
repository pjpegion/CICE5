subroutine fill_vertices(nv,jbeg,jend,ni,iVert,jVert,lat,lon,latvert,lonvert)

  implicit none

                            integer, intent( in) :: nv,nj,ni,jbeg,jend
                            integer, intent( in) :: iVert(nv), jVert(nv)
  real(kind=8), dimension(ni,nj),    intent( in) ::  lat, lon

  real(kind=8), dimension(ni,nj,nv), intent(out) :: latvert,lonvert

  integer :: ii,jj


  do j = jbeg,jend
   do i = 1,ni
    do n = 1,nv
      ii = i + iVert(n); jj = j + jVert(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latvert(i,j,n)   = lat(ii,jj)
      lonvert(i,j,n)   = lon(ii,jj)
    enddo
   enddo
  enddo
end subroutine fill_vertices
 
subroutine fill_bottom((nv,ni,iVert,jVert,lat,lon,latvert,lonvert)

  implicit none

                            integer, intent( in) :: nv,nj,ni
                            integer, intent( in) :: iVert(nv), jVert(nv)
  real(kind=8), dimension(ni,nj),    intent( in) ::  lat, lon

  real(kind=8), dimension(ni,nj,nv), intent(out) :: latvert,lonvert

  ! fill in grid bottom (j=1) 
  ! vertices 1,2 are available
  ! vertices 3,4 must be set manually
     j = 1
   do i = 1,ni
    do n = 1,2
      ii = i + iVert(n); jj = j + jVert(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latvert(i,j,n)   = lat(ii,jj)
      lonvert(i,j,n)   = lon(ii,jj)
    enddo
      lonvert(i,j, 3) = lonvert(i,j,2)
      lonvert(i,j, 4) = lonvert(i,j,1)
      latvert(i,j, 3) = -90.0d0
      latvert(i,j, 4) = -90.0d0
   enddo
  endif
end subroutine fill_bottom

subroutine fill_top(nv,ni,iVert,jVert,lat,lon,latvert,lonvert,xlon,xlat)

  implicit none

                            integer, intent( in) :: nv,nj,ni
                            integer, intent( in) :: iVert(nv), jVert(nv)
  real(kind=8), dimension(ni,nj),    intent( in) ::   lat,  lon
  real(kind=8), dimension(ni,        intent( in) ::  xlat, xlon

  real(kind=8), dimension(ni,nj,nv), intent(out) :: latvert,lonvert

  ! fill in grid top (j=nj)
  ! vertices 3,4 are available
  ! vertices 1,2 must be set manually using 'across seam' values
      j = nj
   do i = 1,ni
    do n = 3,4
      ii = i + iVert(n); jj = j + jVert(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latvert(i,j,n)  = lat(ii,jj)
      lonvert(i,j,n)  = lon(ii,jj)
    enddo
    do n = 1,2
      ii = i + iVert(n)
      if(ii .eq.    0)ii = ni
      if(ii .eq. ni+1)ii = 1
      latvert(i,j,n)  = xlat(ii)
      lonvert(i,j,n)  = xlon(ii)
   enddo
      !latCv_vert(i,j, 1) = latCv_vert(i,j,4)
      !latCv_vert(i,j, 2) = latCv_vert(i,j,3)
      !lonCv_vert(i,j, 1) = lonCv_vert(i,j,4)+240.d0
      !lonCv_vert(i,j, 2) = lonCv_vert(i,j,3)+240.d0
  endif
end subroutine fill_top
