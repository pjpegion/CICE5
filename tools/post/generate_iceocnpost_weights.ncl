;------------------------------------------------------------------
; Denise.Worthen@noaa.gov  (Feb 2019)
;
; This script will generate the weights required for both
; unstaggering of velocities on the MOM6 and CICE tripole grids
; as well as bilinear and conservative regridding to pre-defined
; rectilinear grids. Weights will be generated using NCL's ESMF
; regridding script which relies on ESMF_RegridWeightGen.
;
; Prior to running this script, a file containing all the stagger 
; locations and their vertices must be created. This file
; should be created by reading the MOM6 supergrid. This script
; assumes the file exists and is named "tripole.mx025.nc". 
; 
; Generating accurate conservative weights requires that the 
; vertices of the center T grid cell be given. These values
; are available in the tripole.mx025.nc file. 
;
; The tripole.mx025.nc file is created using the gen_fixgrid.F90
; and associated code. 
;
; Unless the underlying tripole grid changes, this script needs to be
; run one time only.
;
; If the environment variable NumProc is set prior to running this
; script, the weight generation will use NumProc to create the weights.
;
; This file still has to be loaded manually
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
; substitute a curvilinear_to_SCRIP.ncl which does not require
; re-ordering the grid
load "my_curvilinear_to_SCRIP.ncl"
;------------------------------------------------------------------
begin

    ; the weights need to generated one time only. The following
    ; flags will re-generate the chosen weights file
    generate_unstagger_weights = True
    generate_dstgrid_SCRIP_files = True
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

     Opt@GridCornerLon = gf->$cornerlonname$
     Opt@GridCornerLat = gf->$cornerlatname$
     print("using locations "+srclatname+"  "+\
                              srclonname+"  "+\
                           cornerlatname+"  "+\
                           cornerlonname)

    ; substitue local version of curvilinear_to_SCRIP here
    my_curvilinear_to_SCRIP(srcfile, gf->$srclatname$, gf->$srclonname$, Opt)
    delete(Opt) 
   end do 

   ; close
   delete(gf)

;----------------------------------------------------------------------
; generate the weights used to unstagger the grids
;----------------------------------------------------------------------

     dstfile = nemsrc+"Ct_SCRIP.nc"
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

     srcfile = nemsrc+"Ct_SCRIP.nc"
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

