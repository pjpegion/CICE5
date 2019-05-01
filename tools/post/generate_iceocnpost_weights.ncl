;------------------------------------------------------------------
; Denise.Worthen@noaa.gov (Feb 2019)
;
; Post-processing of MOM6/CICE5 tripole grid output
; to regular lat/lon grids 
;
; Prior to running this script, a file containing all the stagger 
; locations and their vertices must be created. This file
; should be created by reading the MOM6 supergrid. This script
; assumes the file exists and is named "tripole.mx025.nc". 
;
; Note: the descriptive text below assumes fortran type indexing
; where the variables are indexed as (i,j) and indices start at 1
; NCL indices are (j,i) and start at 0
;
; The post involves these steps
;
;  a) unstaggering velocity points
;     MOM6 is on an Arakawa C grid. MOM6 refers to these 
;     locations as "Ct" for the centers and "Cu", "Cv" 
;     "Bu" for the left-right, north-south and corner 
;     points, respectively.
;
;     The indexing scheme in MOM6 is as follows:  
;
;        Cv@i,j
;     ----X------X Bu@i,j
;                |
;                |
;       Ct@i,j    |
;         X      X Cu@i,j
;                |
;                |
;                |
;
;     CICE5 is on an Arakawa B grid. CICE5 refers to these
;     locations as TLAT,TLON for the centers and ULAT,ULON
;     for the corners
;  
;     In UFS, the CICE5 grid has been created using the MOM6
;     supergrid file. Therefore, all grid points are consistent
;     between the two models. 
;   
;     In the following, MOM6's nomenclature will be followed,
;     so that CICE5's U-grid will be referred to as "Bu".
;
;  b) rotation of tripole vectors to East-West
;     MOM6 and CICE6 both output velocties on their native
;     velocity points. For MOM6, that is u-velocities on the
;     Cu grid and v-velocites on the Cv grid. For CICE5, it is
;     both u and v-velocities on the Bu grid.
;   
;     The rotation angle for both models are defined at center
;     grid points; therefore the velocities need to be first 
;     unstaggered before rotation. MOM6 and CICE5 also define 
;     opposite directions for the rotations. Finally, while the
;     grid points are identical between the two models, CICE5 
;     calculates the rotation angle at center grid points by 
;     averaging the four surrounding B grid points. MOM6 derives
;     the rotation angle at the center directly from the latitude
;     and longitude of the center grid points. The angles are therefor
;     not identical between the two grids.
; 
; c) conservative regridding of some fields
;    Fields such as ice concentration which are inherently an area
;    weighted value require conservative regridding. Most other
;    variables are state variables and can be regridded using
;    bilinear weighting.
;
;    Generating accurate conservative weights requires that the 
;    vertices of the center T grid cell be given. These values
;    are available in the tripole.mx025.nc file. 
;
; An efficient way to accomplish the unstaggering of velocities 
; is to determine the bilinear interpolation weights between grid 
; points of the Arakawa C grid and the center grid points. For
; example, a set of weights can be generated one time for Cu->T
; and reused thereafter.
;
; An efficient way to regrid from the tripole to rectilinear
; grids is to generate interpolation weights for both bilinear  
; and conservative regridding. These weights can also be generated 
; once and reused thereafter. Bilinear weights generated for the first 
; vertical level can be used on other levels (where the masking 
; changes) by utilizing the correct masking procedure.
; For this purpose, it is useful to generate a SCRIP file for the
; Ct grid which does not contain a mask.
;
; This script will generate the weights required for both
; unstaggering of velocities on the MOM6 and CICE tripole grids
; as well as bilinear and conservative regridding to pre-defined
; regular rectilinear grids. Weights will be generated using NCL's ESMF
; regridding tools.
;
; This file still has to be loaded manually
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
; substitute a curvilinear_to_SCRIP.ncl which does not require
; re-ordering the grid
load "my_curvilinear_to_SCRIP.ncl"
;------------------------------------------------------------------
begin

    ; the weights need to generated one time only. The following
    ; flags will re-generate the chose weights file
    generate_unstagger_weights = False
    generate_dstgrid_SCRIP_files = False
    generate_tripole_rect_weights = True

    ; specify a location to use
        nemsrc     = "/scratch4/NCEPDEV/ocean/save/Denise.Worthen/NEMS_INPUT0.1/ocnicepost/"
    ; interpolation methods
    methods        = (/"bilinear" ,"conserve"/)

    ; source tripole grid  
            srcloc = nemsrc
          gridfile = "tripole.mx025.nc"
       staggerlocs = (/"Ct", "Cu", "Cv", "Bu"/)

    ; destintion grids 
       dstgrds     = (/    "0p25",     "0p5",  "1p0"/)
       dstres      = (/ "0.25deg",  "0.5deg", "1deg"/)

       dstlonbeg   = (/  0.125,    0.25,   0.50/)
       dstlonend   = (/359.875,  359.75, 359.50/)
       dstlatbeg   = (/-89.875,  -89.75, -89.50/)
       dstlatend   = (/ 89.875,   89.75,  89.50/)

        dstxdim    = (/1440,    720,   360/)
        dstydim    = (/ 720,    360,   180/)

   if(generate_unstagger_weights)then
;----------------------------------------------------------------------
; for each stagger grid, first generate SCRIP grid files which will be
; used in the gen_weights call
;----------------------------------------------------------------------

      gf = addfile(nemsrc+gridfile,"r")
   do ii = 0,dimsizes(staggerlocs)-1

    Opt                   = True
    Opt@ForceOverwrite    = True
    Opt@PrintTimings      = True
    Opt@InterpMethod      = "bilinear"
    Opt@Debug             = True
    Opt@Check             = True

    print("working on stagger loc "+staggerlocs(ii))
          srcfile = nemsrc+staggerlocs(ii)+"_SCRIP.nc"
       srclonname = "lon"+staggerlocs(ii)
       srclatname = "lat"+staggerlocs(ii)
    cornerlonname = "lon"+staggerlocs(ii)+"_vert"
    cornerlatname = "lat"+staggerlocs(ii)+"_vert"

    if(staggerlocs(ii) .eq. "Ct")then
     Opt@GridMask      = gf->wet
     print("using grid mask wet")
    end if
     Opt@GridCornerLon = gf->$cornerlonname$
     Opt@GridCornerLat = gf->$cornerlatname$
     print("using locations "+srclatname+"  "+\
                              srclonname+"  "+\
                           cornerlatname+"  "+\
                           cornerlonname)

    my_curvilinear_to_SCRIP(srcfile, gf->$srclatname$, gf->$srclonname$, Opt)
    delete(Opt) 
   end do 

    ; generate a Ct SCRIP file with no masking, this is for generating the
    ; interpolation mask at different vertical ocean model levels
    Opt                   = True
    Opt@ForceOverwrite    = True
    Opt@PrintTimings      = True
    Opt@InterpMethod      = "bilinear"
    Opt@Debug             = True
    Opt@Check             = True

    ii = 0
    print("working on stagger loc "+staggerlocs(ii))
          srcfile = nemsrc+staggerlocs(ii)+"_SCRIP_nomask.nc"
       srclonname = "lon"+staggerlocs(ii)
       srclatname = "lat"+staggerlocs(ii)
    cornerlonname = "lon"+staggerlocs(ii)+"_vert"
    cornerlatname = "lat"+staggerlocs(ii)+"_vert"

     Opt@GridCornerLon = gf->$cornerlonname$
     Opt@GridCornerLat = gf->$cornerlatname$
     print("using locations "+srclatname+"  "+\
                              srclonname+"  "+\
                           cornerlatname+"  "+\
                           cornerlonname)

    my_curvilinear_to_SCRIP(srcfile, gf->$srclatname$, gf->$srclonname$, Opt)
    delete(Opt)

   ; close
   delete(gf)

;----------------------------------------------------------------------
; generate the weights used to unstagger the grids
;----------------------------------------------------------------------

     dstfile = nemsrc+"Ct_SCRIP_nomask.nc"
    do ii = 1,dimsizes(staggerlocs)-1

     Opt                   = True
     Opt@ForceOverwrite    = True
     Opt@PrintTimings      = True
     Opt@InterpMethod      = "bilinear"
     Opt@Debug             = True
     Opt@Check             = True

     srcfile = nemsrc+staggerlocs(ii)+"_SCRIP.nc"
     ; this is the file containing the weights
     wgtfile = nemsrc+"tripole.mx025."+staggerlocs(ii)+".to.Ct.bilinear.nc"

     ESMF_regrid_gen_weights(srcfile,dstfile,wgtfile,Opt)
     delete(Opt)
    end do
    ;close
    delete(dstfile)
    delete(srcfile)
    delete(wgtfile)
   end if

   if(generate_dstgrid_SCRIP_files)then
;----------------------------------------------------------------------
;
;----------------------------------------------------------------------

     do ii = 0,dimsizes(dstgrds)-1

     ; rectilinear destination
     Opt                   = True
     Opt@ForceOverwrite    = True
     Opt@PrintTimings      = True
     Opt@Debug             = True
     Opt@Check             = True

     lons = fspan(dstlonbeg(ii),dstlonend(ii),dstxdim(ii))
     lats = fspan(dstlatbeg(ii),dstlatend(ii),dstydim(ii))
     ;print(dimsizes(lons))
     ;print(dimsizes(lats))
     Opt@LLCorner          = (/ -90.d,   0.d/)
     Opt@URCorner          = (/  90.d, 360.d/)

     Opt@DstGridLat = lats
     Opt@DstGridLon = lons
  
     dstfile = nemsrc+"rect."+dstgrds(ii)+"_SCRIP.nc"
     print("generating dstfile "+dstfile)
     rectilinear_to_SCRIP(dstfile,lats,lons,Opt)
     delete(lons)
     delete(lats)
     delete(Opt)
    end do
   end if

   if(generate_tripole_rect_weights)then
;----------------------------------------------------------------------
; generate SCRIP file for each rectilinear grid and then the weight file
; for tripole Ct->rect
;----------------------------------------------------------------------

    srcfile = nemsrc+"Ct_SCRIP.nc"
    do jj = 0,dimsizes(methods)-1
      Opt                   = True
      Opt@ForceOverwrite    = True
      Opt@PrintTimings      = True
      Opt@InterpMethod      = methods(jj)
      Opt@Debug             = True
      Opt@Check             = True

     srcfile = nemsrc+"Ct_SCRIP_nomask.nc"
     do ii = 0,dimsizes(dstgrds)-1
      dstfile = nemsrc+"rect."+dstgrds(ii)+"_SCRIP.nc"
      ;this file contains the weights
      wgtfile = nemsrc+"tripole.mx025.Ct.to.rect."+dstgrds(ii)+"."+methods(jj)+".nc"
      ESMF_regrid_gen_weights(srcfile,dstfile,wgtfile,Opt)
     end do

    end do

   end if
exit
end

