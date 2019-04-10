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
;    are available in the tripole.grid.mx025.nc file. 
;
; An efficient way to accomplish the unstaggering of velocities 
; is to determine the bilinear interpolation weights between grid 
; points of the Arakawa C grid and the center grid points. For
; example, a set of weights can be generated one time for Cu->T
; and reused thereafter.
;
; An efficient way to regrid from the tripole to regular
; lat/lon grids is to generate interpolation weights for both bilinear  
; and conservative regridding. These weights can also be generated 
; once and reused thereafter. Bilinear weights generated for the first 
; vertical level can be used on other levels (where the masking 
; changes) by utilizing the correct masking procedure.
;
; This script will generate the weights required for both
; unstaggering of velocities on the MOM6 and CICE tripole grids
; as well as bilinear and conservative regridding to pre-defined
; regular latlon grids. Weights will be generated using NCL's ESMF
; regridding tools.
;
; This file still has to be loaded manually
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
;------------------------------------------------------------------
begin

    generate_unstagger_weights = False
    generate_tripole_latlon_weights = True

;---Specify a location to use
        nemsrc     = "/scratch4/NCEPDEV/ocean/save/Denise.Worthen/NEMS_INPUT0.1/ocnicepost/"
;---Interpolation methods
    methods        = (/"bilinear" ,"conserve"/)
;---Source tripole grid location 
            srcloc = nemsrc
          gridfile = "tripole.mx025.nc"
       staggerlocs = (/"Ct", "Cu", "Cv", "Bu"/)
 
   if(generate_unstagger_weights)then
;----------------------------------------------------------------------
; for each stagger grid, first generate SCRIP grid files which will be
; used in the gen_weights call
;----------------------------------------------------------------------

    Opt                   = True
    Opt@ForceOverwrite    = True
    Opt@PrintTimings      = True
    Opt@Debug             = True
    Opt@InterpMethod      = "bilinear"

      gf = addfile(nemsrc+gridfile,"r")
   do ii = 0,dimsizes(staggerlocs)-1
    print("working on stagger loc "+staggerlocs(ii))
          srcfile = nemsrc+staggerlocs(ii)+"_SCRIP.nc"
       srclonname = "lon"+staggerlocs(ii)
       srclatname = "lat"+staggerlocs(ii)
    cornerlonname = "lon"+staggerlocs(ii)+"_vert"
    cornerlatname = "lat"+staggerlocs(ii)+"_vert"

    Opt@GridCornerLon = gf->$cornerlonname$
    Opt@GridCornerLat = gf->$cornerlatname$

    curvilinear_to_SCRIP(srcfile, gf->$srclatname$, gf->$srclonname$, Opt) 
   end do 
   ; close
   delete(gf)

;----------------------------------------------------------------------
; Generate the weights for unstaggering the grids; these are all
; bilinear
;----------------------------------------------------------------------

     dstfile = nemsrc+"Ct_SCRIP.nc"
    do ii = 1,dimsizes(staggerlocs)-1
     srcfile = nemsrc+staggerlocs(ii)+"_SCRIP.nc"
     ; this is the file containing the weights
     wgtfile = nemsrc+"tripole.mx025."+staggerlocs(ii)+".to.Ct.bilinear.nc"

     ESMF_regrid_gen_weights(srcfile,dstfile,wgtfile,Opt)
    end do
    delete(Opt)
    ;close
    delete(dstfile)
    delete(srcfile)
    delete(wgtfile)
   end if

;----------------------------------------------------------------------
; Generate the weights to regrid from tripole Ct to latlon using
; bilinear and conservative regridding
;----------------------------------------------------------------------

;---Destination grid and weights location
            dstloc = nemsrc
;---Destination grid sizes and name
       dstgrds     = (/    "0p25",     "0p5",  "1p0"/)
       dstres      = (/ "0.25deg",  "0.5deg", "1deg"/)

   if(generate_tripole_latlon_weights)then
;----------------------------------------------------------------------
; generate SCRIP file for each latlon grid and then the weight file
; for tripole Ct->latlon
;----------------------------------------------------------------------

    Opt                   = True
    Opt@ForceOverwrite    = True
    Opt@PrintTimings      = True
    Opt@Debug             = True

    ;ii = 2
    do ii = 0,dimsizes(dstgrds)-1
     Opt@LLCorner       = (/ -90.d,   0.d/)
     Opt@URCorner       = (/  90.d, 360.d/)
    
     dstfile = nemsrc+"latlon."+dstgrds(ii)+"_SCRIP.nc"
     latlon_to_SCRIP(dstfile,dstres(ii),Opt)
     
     jj=1
     srcfile = nemsrc+"Ct_SCRIP.nc"
     ;do jj = 0,dimsizes(methods)-1
       Opt@InterpMethod  = methods(jj)
      ; this is the file containing the weights
      wgtfile = nemsrc+"tripole.mx025.Ct.to.latlon."+dstgrds(ii)+"."+methods(jj)+".nc"
      ESMF_regrid_gen_weights(srcfile,dstfile,wgtfile,Opt)
      delete(Opt@InterpMethod)
      delete(Opt@GridMask)
     ;end do
    end do

   end if
exit
end

