;************************************************

  load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
  load "unstaggerCICE5.ncl"

;************************************************
begin

;************************************************
; off-line generated weights for regridding
; this should be the same regrid weights file 
; from MOM6. The T-grid of CICE is created from 
; the MOM6 supergrid file so the T points are
; identical. However, the angles on MOM6's T-grid
; and CICE's T-grid are are not identical. There
; is a sign convention change between MOM6 and CICE;
; also CICE creates the ANGLET from the 
; area-weighted mean of the surrounding 4 U-grid
; angle values
;************************************************

   dirname="/scratch3/NCEPDEV/stmp2/Denise.Worthen/BenchPostTest/"
   ;res="025"
   res="05"
   if(res .eq. "025")then
    outres="_0p25x0p25_CICE"
   end if
   if(res .eq. "05")then
    outres="_0p5x0p5_CICE"
   end if
   system("/bin/mkdir -p "+dirname+res)

   wgtsrc  = "/scratch4/NCEPDEV/ocean/save/Bhavani.Rajan/post/"
   wgtfile = "remap.grid_t."+res+"deg.bil.nc"
      wgts = wgtsrc + wgtfile

;************************************************
; make a list of the files from the run
;************************************************

   icefilelist = systemfunc("ls "+dirname+"ice*.nc") 
          icef = addfiles(icefilelist,"r")
        nfiles = dimsizes(icefilelist)
    ;print(nfiles)
    
   ; use the first file to retrieve grid variables
    tmask = icef[0]->tmask
    tarea = icef[0]->tarea
    uarea = icef[0]->uarea
  ; in radians, on t-grid
   rotang = icef[0]->ANGLET

     dims = dimsizes(tmask)
     jdim = dims(0)
     idim = dims(1)
   print("tripole grid resolution "+idim+"  "+jdim)

;************************************************
; the variables to be regridded 
;************************************************

   varlist = (/   "hi_h" \
              ,   "hs_h" \
              , "Tsfc_h" \
              , "aice_h" \
              ,  "sst_h" \
              , "uvel_h" \
              , "vvel_h" \
              /)

   nvars = dimsizes(varlist)
   print(nvars+" to be re-gridded")
   ; an array to hold the tripole variables
   ; each file contains a single timestep
   vars = new((/nvars,1,jdim,idim/),float)

   ; the u-vel and v-vel index
   do ii = 0,nvars-1
    varname = varlist(ii)
    if(varname .eq. "uvel_h")then
     uindx = ii
    end if
    if(varname .eq. "vvel_h")then
     vindx = ii
    end if
   end do
   print("u velocity at index "+uindx)
   print("v velocity at index "+vindx)

;************************************************
; loop over files
;************************************************

   ;ii = 0
   ;do ii = 0,nfiles-1
   do ii = 0,1
    ; create output file name
    outfile = icefilelist(ii)
    ; find the index where the actual filename begins
     indstr = str_index_of_substr(outfile, "/", -1)
    ; insert a sub-directory name
    outfile = str_insert(outfile,res+"/",indstr+1)
    ; insert a string indicating regridding resolution
    outfile = str_insert(outfile,outres,-4)
    ; find the index where the actual filename begins
     indstr = str_index_of_substr(outfile, "ice", -1)
    ; insert an 'r' after the leading occurence of 'ice' 
    ; to indicate regular grid
    outfile = str_insert(outfile,"r",indstr+3)
    system("/bin/rm -f "+outfile)

         df = addfile(outfile,"c")
    filedimdef(df,"time",-1,True)   ; make time an UNLIMITED dimension 
    print("regridded files will be written to "+outfile)

      time = icef[ii]->time
    ;retrieve the tripole variabes
    do ll = 0,nvars-1
            varname = varlist(ll)
     print(varname)
     vars(ll,:,:,:) = icef[ii]->$varname$(:,:,:)
    end do

    ; unstagger and rotate the velocities
    
    ;printVarSummary(aice_h)
    ; delete the attributes which are no longer valid after regridding
    ;delete_VarAtts(aice_h,(/"coordinates", "cell_measures"/))
    ;delete_VarAtts(  hi_h,(/"coordinates", "cell_measures"/))
    ;delete_VarAtts(uvel_h,(/"coordinates", "cell_measures"/))
    ;delete_VarAtts(vvel_h,(/"coordinates", "cell_measures"/))
 
    ;2d copy for unstaggering; copies metadata
    utmp = vars(uindx,0,:,:)
    vtmp = vars(vindx,0,:,:)
     ; replace missing values before unstaggering
     utmp = where(ismissing(utmp),0.0,utmp)
     vtmp = where(ismissing(vtmp),0.0,vtmp)
     ; move to T-grid
       ut = unstaggerCICE5(utmp,uarea,tarea)
       vt = unstaggerCICE5(vtmp,uarea,tarea)
     ; copy the remaining attributes
     copy_VarAtts(utmp,ut)
     copy_VarAtts(vtmp,vt)
     ;delete(utmp)
     ;delete(vtmp)

     ; rotate
     utmp = cos(rotang)*ut - sin(rotang)*vt
     vtmp = sin(rotang)*ut + cos(rotang)*vt
     ; change attribute to indicate these are now rotated velocities
     utmp@long_name=str_sub_str(utmp@long_name,"(x)","(zonal)")
     vtmp@long_name=str_sub_str(vtmp@long_name,"(y)","(merid)")
     ; re-mask
     utmp = mask(utmp, tmask .eq. 1.0, utmp@_FillValue)
     vtmp = mask(vtmp, tmask .eq. 1.0, utmp@_FillValue)

     ; re-insert
     vars(uindx,:,:,:) = utmp(:,:)
     vars(vindx,:,:,:) = vtmp(:,:)
     ; copy coords from T-grid 
     copy_VarCoords(vars(0),ars[uindx])
     copy_VarCoords(varlist,vars[vindx])
     ; add time coord
     ;df->time = time
     print("regridding file "+ii)
     ; regrid all variables using weights file
     ;df->hi_h   = ESMF_regrid_with_weights(  hi_h, wgts, False)
     ;df->hs_h   = ESMF_regrid_with_weights(  hs_h, wgts, False)
     ;df->Tsfc_h = ESMF_regrid_with_weights(Tsfc_h, wgts, False)
     ;df->aice_h = ESMF_regrid_with_weights(aice_h, wgts, False)
     ;df->sst_h  = ESMF_regrid_with_weights( sst_h, wgts, False)
    ; df->uvel_h = ESMF_regrid_with_weights(uvel_h, wgts, False)
    ; df->vvel_h = ESMF_regrid_with_weights(vvel_h, wgts, False)
    ; delete(df)

   end do ; nfiles

exit
end
