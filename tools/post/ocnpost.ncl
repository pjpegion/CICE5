;------------------------------------------------------------------
; Denise.Worthen@noaa.gov (Feb 2019)
;
; Post-processing of MOM6/CICE5 tripole grid output
; to regular lat/lon grids 
;
; This script will regrid MOM6 ocean output on the
; tripole grid to destination a set of rectilinear
; grids using pre-computed ESMF weights to regrid the
; listed fields to the destination grid and write the 
; results to a new netCDF file
;
  load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
  load "ocnpost_add_functions.ncl"

;----------------------------------------------------------------------
begin

    output_masks = True

   ; specify a location to use
       nemsrc     = "/scratch4/NCEPDEV/ocean/save/Denise.Worthen/NEMS_INPUT0.1/ocnicepost/"
   ; interpolation methods
   methods        = (/"bilinear" ,"conserve"/)
   ; ocean model output location 
           dirsrc = "/scratch3/NCEPDEV/stmp2/Denise.Worthen/BM1_ocn/"

   ; destination grid sizes and name
      dsttype     = (/"rect."/)
      dstgrds     = (/"0p25", "0p5", "1p0"/)

   ; variables to be regridded with the native tripole stagger location
   ; and dimensionality 

   varlist = (/ (/  "SSH", "Ct", "bilinear", "2"/) \
               ,(/  "SST", "Ct", "bilinear", "2"/) \
               ,(/  "SSS", "Ct", "bilinear", "2"/) \
               ,(/"speed", "Ct", "bilinear", "2"/) \
               ,(/ "temp", "Ct", "bilinear", "3"/) \
               ,(/   "so", "Ct", "bilinear", "3"/) \
             /)
   ;print(varlist)

    ; vectors to be regridded with the native tripole stagger location
    ; and dimensionality 

   nvpairs = 2
   veclist = new( (/nvpairs,4,2/),"string")
   veclist = (/ (/ (/"SSU", "SSV"/), (/"Cu", "Cv"/), (/"bilinear", "bilinear"/), (/"2", "2"/) /) ,\
                (/ (/ "uo",  "vo"/), (/"Cu", "Cv"/), (/"bilinear", "bilinear"/), (/"3", "3"/) /)  \
             /)
   ;print(veclist)

     dims = dimsizes(varlist)
    nvars = dims(0)
   delete(dims)

   begTime = get_cpu_time()
;----------------------------------------------------------------------
; make a list of the directories and files from the run 
;----------------------------------------------------------------------

   idate = "20120101"

   ocnfilelist = systemfunc("ls "+dirsrc+"gfs."+idate+"/00/"+"ocn*.nc") 
          ocnf = addfiles(ocnfilelist,"r")
        nfiles = dimsizes(ocnfilelist)

  ; get the rotation angles and vertical grid from the first file
    sinrot = ocnf[0]->sinrot
    cosrot = ocnf[0]->cosrot
       z_l = ocnf[0]->z_l
       z_i = ocnf[0]->z_i
     nlevs = dimsizes(z_l)

  ; get a 2 and 3 dimensional fields for creating the interpolation masks
  mask2d = where(ismissing(ocnf[0]->SST),  1.0, 0.0)
  mask3d = where(ismissing(ocnf[0]->temp), 1.0, 0.0)
  ;printVarSummary(mask2d)
  ;printVarSummary(mask3d)

;----------------------------------------------------------------------
; loop over the output resolutions 
;----------------------------------------------------------------------

      jj = 2
      ii = 0

   do jj = 0,dimsizes(dstgrds)-1
    outres = "_"+dstgrds(jj)+"x"+dstgrds(jj)+"_MOM6"

   ; regrid a field to obtain the output xy dimensions
    wgtsfile = nemsrc+"tripole.mx025.Ct.to."+dsttype+dstgrds(jj)+".bilinear.nc"
          tt = ESMF_regrid_with_weights(sinrot,wgtsfile,False)
        tt!0 = "lat"
        tt!1 = "lon"
         lat = tt&lat
         lon = tt&lon
        dims = dimsizes(tt)
        nlat = dims(0)
        nlon = dims(1)
     print("fields will be remapped to destination grid size "\
           +nlon+"  "+nlat)
  
     delete(tt)
     delete(dims)

   ; regrid the masks to obtain the interpolation masks. keep  the time
   ; dimension (even though they are constant) 
    wgtsfile = nemsrc+"tripole.mx025.Ct.to."+dsttype+dstgrds(jj)+".bilinear.nc"
    rgmask2d = new((/1,nlat,nlon/),typeof(mask2d),default_fillvalue(typeof(mask2d)))
    rgmask3d = new((/1,nlevs,nlat,nlon/),typeof(mask3d),default_fillvalue(typeof(mask2d)))
  
    ; the mask2d,mask3d contain 1's on land and 0's at valid points.
    ; wherever the regridded mask > 0, land values have crept into the
    ; field. remapped fields can then be masked with the interpolation
    ; mask 
    rgmask2d = ESMF_regrid_with_weights(mask2d, wgtsfile,False)
    rgmask3d = ESMF_regrid_with_weights(mask3d, wgtsfile,False)

    ; create the interpolation mask
    rgmask2d = where(rgmask2d .gt. 0.0, rgmask2d@_FillValue, 1.0)
    rgmask3d = where(rgmask3d .gt. 0.0, rgmask3d@_FillValue, 1.0)

;----------------------------------------------------------------------
; loop over each file in the ocnfilelist
;----------------------------------------------------------------------

   ;do ii = 0,nfiles-1
     infile = ocnfilelist(ii)
    ;outfile = create_outfile_name(infile,outres,dstgrds(jj))
    outfile = "test"+outres+".nc"
       time = ocnf[0]->time
    delete(time@bounds)

;----------------------------------------------------------------------
; set up the output netcdf file
;----------------------------------------------------------------------
 
    system("/bin/rm -f " + outfile)    ; remove if exists
    outcdf  = addfile (outfile, "c")  ; open output file

    ; explicitly declare file definition mode. Improve efficiency.
    setfileoption(outcdf,"DefineMode",True)

    ; create global attributes of the file
    fAtt               = True            ; assign file attributes
    fAtt@creation_date = systemfunc ("date")        
    fileattdef( outcdf, fAtt )           ; copy file attributes    

    ; predefine the coordinate variables and their dimensionality
    dimNames = (/"time", "z_l",   "z_i", "lat", "lon"/)  
    dimSizes = (/ -1   , nlevs, nlevs+1,  nlat,  nlon/) 
    dimUnlim = (/ True , False,   False, False, False/)   
    filedimdef(outcdf,dimNames,dimSizes,dimUnlim)

    ; predefine the the dimensionality of the variables to be written out
    filevardef(outcdf, "time", typeof(time), getvardims(time)) 
    filevardef(outcdf,  "z_l",  typeof(z_l),  getvardims(z_l))                           
    filevardef(outcdf,  "z_i",  typeof(z_i),  getvardims(z_i))                           
    filevardef(outcdf,  "lat",  typeof(lat),  getvardims(lat))                          
    filevardef(outcdf,  "lon",  typeof(lon),  getvardims(lon))                          

    ; Copy attributes associated with each variable to the file
    filevarattdef(outcdf, "time", time)                 
    filevarattdef(outcdf,  "z_l",  z_l)               
    filevarattdef(outcdf,  "z_i",  z_i)               
    filevarattdef(outcdf,  "lat",  lat)             
    filevarattdef(outcdf,  "lon",  lon)            

    do nv = 0,nvars-1
     varname = varlist(nv,0)
     vardims = varlist(nv,3)
     if(vardims .eq. "2")then
      odims = (/"time", "lat", "lon"/)
     else
      odims = (/"time", "z_l", "lat", "lon"/)
     end if
     print("creating variable "+varname+" in file")
     filevardef(outcdf, varname, "float", odims)
     delete(odims)
    end do

    do nv = 0,nvpairs-1
     do nn = 0,1
     vecname = veclist(nv,0,nn)
     vecdims = veclist(nv,3,nn)
     if(vecdims .eq. "2")then
      odims = (/"time", "lat", "lon"/)
     else
      odims = (/"time", "z_l", "lat", "lon"/)
     end if
     print("creating variable "+vecname+" in file")
     filevardef(outcdf, vecname, "float", odims) 
     delete(odims)
     delete(vecdims)
     end do
    end do
 
    if(output_masks)then
     filevardef(outcdf, "mask2d", "float", (/"time", "lat", "lon"/)) 
     filevardef(outcdf, "mask3d", "float", (/"time", "z_l", "lat", "lon"/)) 
    end if

    ; explicitly exit file definition mode.
    setfileoption(outcdf,"DefineMode",False)

    ; write the dimensions to the file
    outcdf->time   = (/time/)     
    outcdf->z_l    = (/z_l/)     
    outcdf->z_i    = (/z_i/)     
    outcdf->lat    = (/lat/)
    outcdf->lon    = (/lon/) 

    if(output_masks)then
     outcdf->mask2d    = (/rgmask2d/) 
     outcdf->mask3d    = (/rgmask3d/) 
    end if

;----------------------------------------------------------------------
;
;----------------------------------------------------------------------

    do nv = 0,nvars-1
     varname = varlist(nv,0)
     vargrid = varlist(nv,1)
     varmeth = varlist(nv,2)
     vardims = varlist(nv,3)
     print(nv+"   "+varname+"  "+vargrid+"  "+varmeth)
    
    ocnvar = ocnf[ii]->$varname$
     ndims = dimsizes(dimsizes(ocnvar))
     ;print(ndims+"   "+dimsizes(ocnvar))

     if(vargrid .ne. "Ct")then
      ; print error if the variable is not on the Ct grid
      print("Variable is not on Ct grid")
      exit
     end if

     ; regrid to dsttype+dstgrd with method
     print("remapping "+varname+" to grid "+dsttype+dstgrds(jj))
     wgtsfile = nemsrc+"tripole.mx025.Ct.to."+dsttype+dstgrds(jj)+"."+varmeth+".nc"

     rgtt = ESMF_regrid_with_weights(ocnvar,wgtsfile,False)
    if(vardims .eq. "2")then
     rgtt = where(ismissing(rgmask2d),ocnvar@_FillValue,rgtt)
    else
     rgtt = where(ismissing(rgmask3d),ocnvar@_FillValue,rgtt)
    end if

     ; enter file definition mode to add variable attributes
     setfileoption(outcdf,"DefineMode",True)
     filevarattdef(outcdf, varname, rgtt)                 
     setfileoption(outcdf,"DefineMode",False)

     outcdf->$varname$   = (/rgtt/)

     delete(ocnvar)
     delete(rgtt)     
    ; nv, loop over number of variables
    end do

;----------------------------------------------------------------------
;
;----------------------------------------------------------------------
   
   do nv = 0,nvpairs-1
     vecnames = veclist(nv,0,:)
     vecgrids = veclist(nv,1,:)
     vecmeth  = veclist(nv,2,:)
     vecdims  = veclist(nv,3,:)
     print(nv+"   "+vecnames+"  "+vecgrids+"  "+vecmeth)

     ; create a vector pair list
     vecpairs = NewList("fifo")
            n = 0
         uvel = ocnf[ii]->$vecnames(n)$
       vecfld = where(ismissing(uvel),0.0,uvel)
        copy_VarAtts(uvel,vecfld)
     print("unstagger "+vecnames(n)+" from "+vecgrids(n)+" to Ct")
     wgtsfile = nemsrc+"tripole.mx025."+vecgrids(n)+".to.Ct.bilinear.nc"
           ut = ESMF_regrid_with_weights(vecfld,wgtsfile,False)

            n = 1
         vvel = ocnf[ii]->$vecnames(n)$
       vecfld = where(ismissing(vvel),0.0,vvel)
        copy_VarAtts(vvel,vecfld)
     print("unstagger "+vecnames(n)+" from "+vecgrids(n)+" to Ct")
     wgtsfile = nemsrc+"tripole.mx025."+vecgrids(n)+".to.Ct.bilinear.nc"
           vt = ESMF_regrid_with_weights(vecfld,wgtsfile,False)

     ListAppend(vecpairs,ut)
     ListAppend(vecpairs,vt)
     ;print(vecpairs)

     ; rotate
     ;vecpairs = rotate_to_EW(vecpairs,cosrot,sinrot,nlevs)
     ;print(vecpairs)

     ; remap
     do n = 0,1
      vecfld = vecpairs[n]
      ; regrid to dsttype+dstgrd with method
      print("remapping "+vecnames(n)+" to grid "+dsttype+dstgrds(jj))
      wgtsfile = nemsrc+"tripole.mx025.Ct.to."+dsttype+dstgrds(jj)+"."+vecmeth(n)+".nc"

       rgtt = ESMF_regrid_with_weights(vecfld,wgtsfile,False)
      if(vecdims(n) .eq. "2")then
       rgtt = where(ismissing(rgmask2d),vecfld@_FillValue,rgtt)
      else
       rgtt = where(ismissing(rgmask3d),vecfld@_FillValue,rgtt)
      end if

      ; enter file definition mode to add variable attributes
      setfileoption(outcdf,"DefineMode",True)
      filevarattdef(outcdf, vecnames(n), rgtt)                 
      setfileoption(outcdf,"DefineMode",False)
    
      outcdf->$vecnames(n)$   = (/rgtt/)
      delete(rgtt)
     end do
      delete([/uvel,vvel,ut,vt,vecfld,vecpairs/])
      delete([/vecnames,vecgrids,vecmeth,vecdims/])
    ; nv, loop over number of vector pairs
    end do

   ; ii, loop over files
   ;end do
  ; jj, loop over destination grids
   delete([/lat,lon,nlon,nlat/])
   delete([/rgmask2d,rgmask3d/])
  end do
  print("One complete ocn file in " + (get_cpu_time() - begTime) + " seconds")
exit
end
