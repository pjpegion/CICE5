module icegriddefs

  implicit none

  integer, parameter :: maxvars = 20

  type icedefs
    character(len=12)   ::  var_name
    character(len=64)   :: long_name
    character(len=12)   :: unit_name
  end type icedefs

  type(icedefs) :: icegrid(maxvars)
  contains

  subroutine ice_typedefine

  integer :: ii = 0

   ii = ii + 1
   icegrid(ii)%var_name  = 'ulon'
   icegrid(ii)%long_name = 'Longitude of corner (Bu) points'
   icegrid(ii)%unit_name = 'radians'

   ii = ii + 1
   icegrid(ii)%var_name  = 'ulat'
   icegrid(ii)%long_name = 'Latitude of corner (Bu) points'
   icegrid(ii)%unit_name = 'radians'

   ii = ii + 1
   icegrid(ii)%var_name  = 'hte'
   icegrid(ii)%long_name = 'Distance between corner (Bu) points, east face'
   icegrid(ii)%unit_name = 'cm'

   ii = ii + 1
   icegrid(ii)%var_name  = 'htn'
   icegrid(ii)%long_name = 'Distance between corner (Bu) points, north face'
   icegrid(ii)%unit_name = 'cm'

   ii = ii + 1
   icegrid(ii)%var_name  = 'angle'
   icegrid(ii)%long_name = 'Angle at corner (Bu) points'
   icegrid(ii)%unit_name = 'radians'

   ii = ii + 1
   icegrid(ii)%var_name  = 'anglet'
   icegrid(ii)%long_name = 'Angle at center (T) points'
   icegrid(ii)%unit_name = 'radians'

   ii = ii + 1
   icegrid(ii)%var_name  = 'lonT'
   icegrid(ii)%long_name = 'Longitude of center (T) points'
   icegrid(ii)%unit_name = 'degrees'

   ii = ii + 1
   icegrid(ii)%var_name  = 'latT'
   icegrid(ii)%long_name = 'Latitude of center (T) points'
   icegrid(ii)%unit_name = 'degrees'

   ii = ii + 1
   icegrid(ii)%var_name  = 'lonCv'
   icegrid(ii)%long_name = 'Longitude of meridional velocity (Cv) points'
   icegrid(ii)%unit_name = 'degrees'

   ii = ii + 1
   icegrid(ii)%var_name  = 'latCv'
   icegrid(ii)%long_name = 'Latitude of meridional velocity (Cv) points'
   icegrid(ii)%unit_name = 'degrees'

   ii = ii + 1
   icegrid(ii)%var_name  = 'lonCu'
   icegrid(ii)%long_name = 'Longitude of zonal velocity (Cu) points'
   icegrid(ii)%unit_name = 'degrees'

   ii = ii + 1
   icegrid(ii)%var_name  = 'latCu'
   icegrid(ii)%long_name = 'Latitude of zonal velocity (Cu) points'
   icegrid(ii)%unit_name = 'degrees'

 end subroutine ice_typedefine
end module icegriddefs
