;************************************************

  load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
  load "unstaggerCICE5.ncl"

;************************************************
begin

;************************************************
; off-line generated weights for regridding
; this is a sample file for testing; you should
; be able to use the same regrid weights file 
; from MOM6. The T-grid of CICE is created from 
; the MOM6 supergrid file so the T points are
; identical. However, the angles on MOM6's T-grid
; and CICE's T-grid are are not identical. There
; is a sign convention change between MOM6 and CICE;
; also CICE creates the ANGLET from the 
; area-weighted mean of the surrounding 4 U-grid
; angle values
;************************************************

   nemsrc  = "/scratch4/NCEPDEV/ocean/save/Denise.Worthen/NEMS_INPUT0.1/"
   wgtsrc  = nemsrc+"regrids/"
   wgtfile = "mom6_rect024.nc"
      wgts = wgtsrc + wgtfile

;************************************************
; the file containing the destination grid used
; to create the above weights file; you'll need
; to supply the specific names in your remap file
;************************************************

        dstgrid = "etopo024_oceanmask.nc"
     dstlatname = "lat"
     dstlonname = "lon"
    dstmaskname = "wet"

      gridfile = addfile(wgtsrc+dstgrid,"r")
          dstlon = gridfile->$dstlonname$
          dstlat = gridfile->$dstlatname$
         dstmask = gridfile->$dstmaskname$

;************************************************
; make a list of the files from the run
; this will 'catch' the initial condition file
; until Bin rewrites the workflow to give the 
; initial condition file a more unique name
;************************************************

   icefilelist = systemfunc("ls "+dirname+"ice*.nc") 
          icef = addfiles(icefilelist,"r")
        nfiles = dimsizes(icefilelist)
   ; print(nfiles)
    
   ; use the first file to retrieve grid variables
    tmask = icef[0]->tmask
    tarea = icef[0]->tarea
    uarea = icef[0]->uarea
  ; in radians, on t-grid
   rotang = icef[0]->ANGLE

;************************************************
; loop over files
;************************************************

   ii = 0
   ;do ii = 0,nfiles-1
    ; create output file name
    outfile = icefilelist(ii)
    ; find the index where the actual filename begins
     indstr = str_index_of_substr(outfile, "/", -1)
    ; insert a sub-directory name
    outfile = str_insert(outfile,"025/",indstr+1)
    ; insert a string indicating regridding resolution
    outfile = str_insert(outfile,"025.",-3)
         df = addfile(outfile,"c")
    filedimdef(df,"time",-1,True)
    print("regridded files will be written to "+outfile)

    aice_h = icef[ii]->aice_h(ii,:,:)
      hi_h = icef[ii]->hi_h(ii,:,:)
    uvel_h = icef[ii]->uvel_h(ii,:,:)
    vvel_h = icef[ii]->vvel_h(ii,:,:)
    ;printVarSummary(aice_h)
    ; delete the attributes which are no longer valid after regridding
    delete_VarAtts(aice_h,(/"coordinates", "cell_measures"/))
    delete_VarAtts(  hi_h,(/"coordinates", "cell_measures"/))
    delete_VarAtts(uvel_h,(/"coordinates", "cell_measures"/))
    delete_VarAtts(vvel_h,(/"coordinates", "cell_measures"/))

     ; replace missing values before unstaggering
     uvel_h = where(ismissing(uvel_h),0.0,uvel_h)
     vvel_h = where(ismissing(vvel_h),0.0,vvel_h)
     ; move to T-grid
       ut = unstaggerCICE5(uvel_h,uarea,tarea)
       vt = unstaggerCICE5(vvel_h,uarea,tarea)
     ; copy the remaining attributes
     copy_VarAtts(uvel_h,ut)
     copy_VarAtts(vvel_h,vt)

     ; rotate
     uvel_h = cos(rotang)*ut - sin(rotang)*vt
     vvel_h = sin(rotang)*ut + cos(rotang)*vt
     ; change attribute to indicate these are now rotated velocities
     uvel_h@long_name=str_sub_str(uvel_h@long_name,"(x)","(zonal)")
     vvel_h@long_name=str_sub_str(vvel_h@long_name,"(y)","(merid)")
  
     ; re-mask
     uvel_h = mask(uvel_h, tmask .eq. 1.0, uvel_h@_FillValue)
     vvel_h = mask(vvel_h, tmask .eq. 1.0, vvel_h@_FillValue)

     ; copy coords from T-grid 
     copy_VarCoords(aice_h,uvel_h)
     copy_VarCoords(aice_h,vvel_h)

     print("regridding file "+ii)
     ; regrid all variables using weights file
     df->aice_h = ESMF_regrid_with_weights(aice_h, wgts, True)
     df->hi_h = ESMF_regrid_with_weights(hi_h, wgts, True)
     df->uvel_h  = ESMF_regrid_with_weights(   uvel_h, wgts, True)
     df->vvel_h  = ESMF_regrid_with_weights(   vvel_h, wgts, True)

   ;end do ; nfiles


exit
end
