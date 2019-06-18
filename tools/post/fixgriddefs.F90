module fixgriddefs

  implicit none

  integer, parameter :: maxvars = 40

  type grid_defs
    character(len=12)   ::  var_name
    character(len=64)   :: long_name
    character(len=20)   :: unit_name
    character(len= 2)   ::  var_type
    character(len=20)   ::  vertices
  end type grid_defs

  type(grid_defs) :: fixgrid(maxvars)
  contains

  subroutine fixgrid_typedefine

  integer :: ii = 0
  
   !default
   fixgrid(:)%var_type  = 'r8'
   fixgrid(:)%vertices  = ' '

   ii = ii + 1
   fixgrid(ii)%var_name  = 'lonCt'
   fixgrid(ii)%long_name = 'Longitude of center (Ct) points'
   fixgrid(ii)%unit_name = 'degrees_east'
   fixgrid(ii)%vertices  = 'lonCt_vert'

   ii = ii + 1
   fixgrid(ii)%var_name  = 'latCt'
   fixgrid(ii)%long_name = 'Latitude of center (Ct) points'
   fixgrid(ii)%unit_name = 'degrees_north'
   fixgrid(ii)%vertices  = 'latCt_vert'

   ii = ii + 1
   fixgrid(ii)%var_name  = 'lonCv'
   fixgrid(ii)%long_name = 'Longitude of meridional velocity (Cv) points'
   fixgrid(ii)%unit_name = 'degrees_east'
   fixgrid(ii)%vertices  = 'lonCv_vert'

   ii = ii + 1
   fixgrid(ii)%var_name  = 'latCv'
   fixgrid(ii)%long_name = 'Latitude of meridional velocity (Cv) points'
   fixgrid(ii)%unit_name = 'degrees_north'
   fixgrid(ii)%vertices  = 'latCv_vert'

   ii = ii + 1
   fixgrid(ii)%var_name  = 'lonCu'
   fixgrid(ii)%long_name = 'Longitude of zonal velocity (Cu) points'
   fixgrid(ii)%unit_name = 'degrees_east'
   fixgrid(ii)%vertices  = 'lonCu_vert'

   ii = ii + 1
   fixgrid(ii)%var_name  = 'latCu'
   fixgrid(ii)%long_name = 'Latitude of zonal velocity (Cu) points'
   fixgrid(ii)%unit_name = 'degrees_north'
   fixgrid(ii)%vertices  = 'latCu_vert'

   ii = ii + 1
   fixgrid(ii)%var_name  = 'lonBu'
   fixgrid(ii)%long_name = 'Longitude of corner (Bu) points'
   fixgrid(ii)%unit_name = 'degrees_east'
   fixgrid(ii)%vertices  = 'lonBu_vert'

   ii = ii + 1
   fixgrid(ii)%var_name  = 'latBu'
   fixgrid(ii)%long_name = 'Latitude of corner (Bu) points'
   fixgrid(ii)%unit_name = 'degrees_north'
   fixgrid(ii)%vertices  = 'latBu_vert'

   ii = ii + 1
   fixgrid(ii)%var_name  = 'lonCt_vert'
   fixgrid(ii)%long_name = 'Longitude Vertices of Ct points'
   fixgrid(ii)%unit_name = 'degrees_east'

   ii = ii + 1
   fixgrid(ii)%var_name  = 'latCt_vert'
   fixgrid(ii)%long_name = 'Latitude Vertices of Ct points'
   fixgrid(ii)%unit_name = 'degrees_north'

   ii = ii + 1
   fixgrid(ii)%var_name  = 'lonCu_vert'
   fixgrid(ii)%long_name = 'Longitude Vertices of Cu points'
   fixgrid(ii)%unit_name = 'degrees_east'

   ii = ii + 1
   fixgrid(ii)%var_name  = 'latCu_vert'
   fixgrid(ii)%long_name = 'Latitude Vertices of Cu points'
   fixgrid(ii)%unit_name = 'degrees_north'

   ii = ii + 1
   fixgrid(ii)%var_name  = 'lonCv_vert'
   fixgrid(ii)%long_name = 'Longitude Vertices of Cv points'
   fixgrid(ii)%unit_name = 'degrees_east'

   ii = ii + 1
   fixgrid(ii)%var_name  = 'latCv_vert'
   fixgrid(ii)%long_name = 'Latitude Vertices of Cv points'
   fixgrid(ii)%unit_name = 'degrees_north'

   ii = ii + 1
   fixgrid(ii)%var_name  = 'lonBu_vert'
   fixgrid(ii)%long_name = 'Longitude Vertices of Bu points'
   fixgrid(ii)%unit_name = 'degrees_east'

   ii = ii + 1
   fixgrid(ii)%var_name  = 'latBu_vert'
   fixgrid(ii)%long_name = 'Latitude Vertices of Bu points'
   fixgrid(ii)%unit_name = 'degrees_north'

 end subroutine fixgrid_typedefine
end module fixgriddefs
