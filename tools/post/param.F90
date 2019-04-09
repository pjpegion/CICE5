module param

  implicit none

  integer, parameter :: ni = 1440, nj = 1080, nv = 4

  ! super-grid source variables
  integer, parameter :: nx  = ni*2, ny  = nj*2

  ! ij offsets moving counter-clockwise around each Ct(i,j)
  integer, parameter, dimension(nv) :: iVertCt = (/0, -1, -1,  0/)
  integer, parameter, dimension(nv) :: jVertCt = (/0,  0, -1, -1/)

  integer, parameter :: ncoord = 2*4             ! 4sets of lat/lon pairs
  integer, parameter :: nverts = 2*4             ! 4sets of lat/lon pairs vertices
  integer, parameter ::  nvars = ncoord + nverts

end module param
