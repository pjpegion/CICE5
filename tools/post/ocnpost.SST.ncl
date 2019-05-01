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
   
   varlist = (/ (/ "sst", "Ct", "bilinear", "2"/) \
             /)
   print(varlist)

     dims = dimsizes(varlist)
    nvars = dims(0)
   delete(dims)

   begTime = get_cpu_time()
;----------------------------------------------------------------------
; make a list of the directories and files from the run 
;----------------------------------------------------------------------

   idate = "20120101"

   ocnfilelist = systemfunc("ls "+dirsrc+"gfs."+idate+"/00/"+"SST*.nc") 
          ocnf = addfiles(ocnfilelist,"r")
        nfiles = dimsizes(ocnfilelist)

  ; get a 2 dimensional fields for creating the interpolation masks
  mask2d = where(ismissing(ocnf[0]->sst),  1.0, 0.0)

;----------------------------------------------------------------------
; loop over the output resolutions 
;----------------------------------------------------------------------

      jj = 1
      ii = 0

   do jj = 0,dimsizes(dstgrds)-1
    outres = "_"+dstgrds(jj)+"x"+dstgrds(jj)+"_MOM6"

   ; regrid a field to obtain the output xy dimensions
    wgtsfile = nemsrc+"tripole.mx025.Ct.to."+dsttype+dstgrds(jj)+".bilinear.nc"
          tt = ESMF_regrid_with_weights(mask2d(0,:,:),wgtsfile,False)
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
  
    ; the mask2d,mask3d contain 1's on land and 0's at valid points.
    ; wherever the regridded mask > 0, land values have crept into the
    ; field. remapped fields can then be masked with the interpolation
    ; mask 
    rgmask2d = ESMF_regrid_with_weights(mask2d, wgtsfile,False)

    ; create the interpolation mask
    rgmask2d = where(rgmask2d .gt. 0.0, rgmask2d@_FillValue, 1.0)

;----------------------------------------------------------------------
; loop over each file in the ocnfilelist
;----------------------------------------------------------------------

   ;do ii = 0,nfiles-1
     infile = ocnfilelist(ii)
    ;outfile = create_outfile_name(infile,outres,dstgrds(jj))
    outfile = "SST"+outres+".nc"
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
    dimNames = (/"time", "lat", "lon"/)  
    dimSizes = (/ -1   ,  nlat,  nlon/) 
    dimUnlim = (/ True , False, False/)   
    filedimdef(outcdf,dimNames,dimSizes,dimUnlim)

    ; predefine the the dimensionality of the variables to be written out
    filevardef(outcdf, "time", typeof(time), getvardims(time)) 
    filevardef(outcdf,  "lat",  typeof(lat),  getvardims(lat))                          
    filevardef(outcdf,  "lon",  typeof(lon),  getvardims(lon))                          

    ; Copy attributes associated with each variable to the file
    filevarattdef(outcdf, "time", time)                 
    filevarattdef(outcdf,  "lat",  lat)             
    filevarattdef(outcdf,  "lon",  lon)            

    ; there is only 1 variable 
     varname = varlist(0)
     vardims = varlist(3)
      odims = (/"time", "lat", "lon"/)
     print("creating variable "+varname+" in file")
     filevardef(outcdf, varname, "float", odims)
     delete(odims)

    ; explicitly exit file definition mode.
    setfileoption(outcdf,"DefineMode",False)

    ; write the dimensions to the file
    outcdf->time   = (/time/)     
    outcdf->lat    = (/lat/)
    outcdf->lon    = (/lon/) 

;----------------------------------------------------------------------
;
;----------------------------------------------------------------------

     varname = varlist(0)
     vargrid = varlist(1)
     varmeth = varlist(2)
     vardims = varlist(3)
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
     rgtt = where(ismissing(rgmask2d),ocnvar@_FillValue,rgtt)

     ; enter file definition mode to add variable attributes
     setfileoption(outcdf,"DefineMode",True)
     filevarattdef(outcdf, varname, rgtt)                 
     setfileoption(outcdf,"DefineMode",False)

     outcdf->$varname$   = (/rgtt/)

     delete(ocnvar)
     delete(rgtt) 
    
   ; ii, loop over files
   ;end do
   ;jj, loop over destination grids
   delete([/lat,lon,nlon,nlat/])
   delete([/rgmask2d/])
  end do
  print("One complete ocn file in " + (get_cpu_time() - begTime) + " seconds")
exit
end
