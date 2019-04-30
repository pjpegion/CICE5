module grdvar

  use param
 
  implicit none

  ! pole locations
  integer(kind=4) :: ipole(2)

  ! land mask
  integer(kind=4), dimension(ni,nj) :: wet

  ! grid stagger locations
  real(kind=8), dimension(ni,nj) :: latCt, lonCt ! lat and lon of T on C-grid
  real(kind=8), dimension(ni,nj) :: latCv, lonCv ! lat and lon of V on C-grid
  real(kind=8), dimension(ni,nj) :: latCu, lonCu ! lat and lon of U on C-grid
  real(kind=8), dimension(ni,nj) :: latBu, lonBu ! lat and lon of corners on C-grid

  ! areas of Ct grid cell
  real(kind=8), dimension(ni,nj) :: areaCt

  ! vertices of each stagger location
  real(kind=8), dimension(ni,nj,nv) :: latCt_vert, lonCt_vert
  real(kind=8), dimension(ni,nj,nv) :: latCu_vert, lonCu_vert
  real(kind=8), dimension(ni,nj,nv) :: latCv_vert, lonCv_vert
  real(kind=8), dimension(ni,nj,nv) :: latBu_vert, lonBu_vert

  integer, dimension(nv) :: iVertBu, iVertCu, iVertCv
  integer, dimension(nv) :: jVertBu, jVertCu, jVertCv

  ! need across seam values of Ct,Cu points to retrieve vertices of Bu and Cv grids
  real(kind=8), dimension(ni) :: xlonCt, xlatCt
  real(kind=8), dimension(ni) :: xlonCu, xlatCu
  ! latitude spacing at bottom of grid
  real(kind=8), dimension(ni) :: dlatBu, dlatCv

end module grdvar
