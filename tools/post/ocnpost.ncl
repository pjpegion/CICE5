;------------------------------------------------------------------
; Denise.Worthen@noaa.gov (Feb 2019)
;
; Post-processing of MOM6/CICE5 tripole grid output
; to regular lat/lon grids 
;
; This script will regrid MOM6 ocean output on the
; tripole grid to destination regular latlon grids
; using pre-computed ESMF weights to unstagger the
; velocities, regrid to the regular grid and then
; write the fields to a netcdf file
;
  load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
  load "ocnpost_add_functions.ncl"

;----------------------------------------------------------------------
begin

;---Specify a location to use
        nemsrc     = "/scratch4/NCEPDEV/ocean/save/Denise.Worthen/NEMS_INPUT0.1/ocnicepost/"
;---Interpolation methods
    methods        = (/"bilinear" ,"conserve"/)
;---Ocean model output location 
            dirsrc = "/scratch3/NCEPDEV/stmp2/Denise.Worthen/BM1_ocn/"

;---Destination grid sizes and name
       dstgrds     = (/    "0p25",     "0p5",  "1p0"/)
       dstres      = (/ "0.25deg",  "0.5deg", "1deg"/)
       ; dstxdim    = (/      1440,       720,   360/)
       ; dstydim    = (/       720,       360,   180/)

;---The variables to be regridded with the native tripole stagger location
;---and dimensionality 

   varlist = (/ (/ "SSH", "Ct", "bilinear", "2"/) \
               ,(/ "SST", "Ct", "bilinear", "2"/) \
               ,(/ "SSS", "Ct", "bilinear", "2"/) \
               ,(/"temp", "Ct", "bilinear", "3"/) \
               ,(/  "so", "Ct", "bilinear", "3"/) \
             /)
   print(varlist)

   nvpairs = 2
   veclist = new( (/nvpairs,4,2/),"string")
   veclist = (/ (/ (/"SSU", "SSV"/), (/"Cu", "Cv"/), (/"bilinear", "bilinear"/), (/"2", "2"/) /) ,\
                (/ (/ "uo",  "vo"/), (/"Cu", "Cv"/), (/"bilinear", "bilinear"/), (/"3", "3"/) /)  \
             /)
   print(veclist)

     dims = dimsizes(varlist)
    nvars = dims(0)
   delete(dims)

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

;----------------------------------------------------------------------
; loop over the output resolutions for each file 
;----------------------------------------------------------------------

      jj = 0
      ii = 0

   ;do jj = 0,dimsizes(dstgrds)-1
    outres = "_"+dstgrds(jj)+"x"+dstgrds(jj)+"_MOM6"
   ;do ii = 0,nfiles-1
     infile = ocnfilelist(ii)
    ;outfile = create_outfile_name(infile,outres,dstgrds(jj))
    outfile = "test.nc"
       time = ocnf[0]->time
    delete(time@bounds)

   ; regrid a field to obtain the output xy dimensions
    wgtsfile = nemsrc+"tripole.mx025.Ct.to.latlon.0p25.bilinear.nc"
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
     filevardef(outcdf, vecname, "float", odims) 
     delete(odims)
     end do
    end do
 
    ; explicitly exit file definition mode. **NOT REQUIRED**
    setfileoption(outcdf,"DefineMode",False)

    outcdf->time   = (/time/)     
    outcdf->z_l    = (/z_l/)     
    outcdf->z_i    = (/z_i/)     
    outcdf->lat    = (/lat/)
    outcdf->lon    = (/lon/) 

    delete(outcdf)
;----------------------------------------------------------------------
;
;----------------------------------------------------------------------

if(1 .eq. 0)then
;----------------------------------------------------------------------
;
;----------------------------------------------------------------------
 
    nv = 0
    ;do nv = 0,nvars-1
     varname = varlist(nv,0)
     vargrid = varlist(nv,1)
     varmeth = varlist(nv,2)
     print(nv+"   "+varname+"  "+vargrid+"  "+varmeth)
    
    ocnvar = ocnf[ii]->$varname$
     ndims = dimsizes(dimsizes(ocnvar))
     print(ndims+"   "+dimsizes(ocnvar))

     wgtsfile = nemsrc+"tripole.mx025."+vargrid+".to.Ct.bilinear.nc"
     ; unstagger the field if required
     tt = where(ismissing(ocnvar),0.0,ocnvar)
     if(vargrid .ne. "Ct")then
        tt = ESMF_regrid_with_weights(tt,wgtsfile,False)
     end if

     ; regrid to latlon
     wgtsfile = nemsrc+"tripole.mx025.Ct.to.latlon.0p25.bilinear.nc"
     ; first get a field of the right dimensionality to use as a mask 
     ; interpolate to latlon; where the resulting field is >0 is where missing values
     ; were interpolated from land into field; mask these out
     if(ndims .eq. 3)then
       tmask = where(ismissing(ocnf[ii]->SSS),1.0,0.0)
      rgmask = ESMF_regrid_with_weights(tmask,wgtsfile,False)
      rtmask = where(rgmask .gt. 0.0,tt@_FillValue,1.0)

       rgtt = ESMF_regrid_with_weights(tt,wgtsfile,False)
       rgtt = where(ismissing(rtmask),tt@_FillValue,rgtt)
     else
      do kk = 0,nlevs-1
         tmask = where(ismissing(ocnf[ii]->so(:,kk,:,:)), 1.0,0.0)
        rgmask = ESMF_regrid_with_weights(tmask,wgtsfile,False)
        rtmask = where(rgmask .gt. 0.0,tt@_FillValue,1.0)

       ttmp = tt(:,kk,:,:)
       rgtt = ESMF_regrid_with_weights(ttmp,wgtsfile,False)
       rgtt = where(ismissing(rtmask),tt@_FillValue,rgtt)
      end do
     end if

     ;if(ndims .eq. 3)then
     ;  tmask = where(ismissing(ocnf[ii]->SSS),1.0,0.0)
     ;   rgtt = regrid_and_mask(tmask,tt,wgtsfile)
     ;else
     ; do kk = 0,nlevs-1
     ;    tmask = where(ismissing(ocnf[ii]->so(:,kk,:,:)), 1.0,0.0)
     ;               ttmp = tt(:,kk,:,:)
     ;     rgtt(:,kk,:,:) = regrid_and_mask(tmask,ttmp,wgtsfile)
     ; end do
     ;end if
    ;printVarSummary(rgtt)
    ; nv
    ;end do
   
    df = addfile("test.nc","c")
    filedimdef(df,"time",-1,True)   ; make time an UNLIMITED dimension
    time = ocnf[ii]->time
    df->sss = rgtt
end if   
 
if(1.eq.0)then
;----------------------------------------------------------------------
;
;----------------------------------------------------------------------
   
    nv = 0
   ;do nv = 0,nvpairs-1
     vecnames = veclist(nv,0,:)
     vecgrids = veclist(nv,1,:)
     print(nv+"   "+vecnames+"  "+vecgrids)

      uvar = ocnf[ii]->$vecnames(0)$
      vvar = ocnf[ii]->$vecnames(1)$
     ndims = dimsizes(dimsizes(uvar))

     ; get a t-grid variable of the right dimensionality to act as regridding mask
     if(ndims .eq. 2)then
      uvtmask = where(ismissing(ocnf[ii]->SSS(0,:,:)),1.0,0.0)
     else
      uvtmask = where(ismissing(ocnf[ii]->so(0,:,:,:)), 1.0,0.0)
     end if

      ;copy metadata and fill missing values prior to unstaggering
      ut = where(ismissing(uvar),0.0,uvar)
      vt = where(ismissing(vvar),0.0,vvar)
      ;print(min(ut)+"   "+max(ut))
      ;print(min(vt)+"   "+max(vt))

     if(ndims .eq. 3)then
       ut = ESMF_regrid_with_weights(ut,nemsrc+"tripole.mx025."+vecgrids(0)+".to.Ct.bilinear.nc",False)
       vt = ESMF_regrid_with_weights(vt,nemsrc+"tripole.mx025."+vecgrids(1)+".to.Ct.bilinear.nc",False)
     ; now rotate; copy metadata first
     urot = ut
     vrot = vt
     ;urot = 
     ;vrot = 
     else
      nlevs = dims(1)
      do kk = 0,nlevs-1
          utmp(:,:) = ut(0,kk,:,:)
          vtmp(:,:) = vt(0,kk,:,:)
       ut(0,kk,:,:) =  ESMF_regrid_with_weights(utmp,nemsrc+"tripole.mx025."+vecgrids(0)+".to.Ct.bilinear.nc",False)
       vt(0,kk,:,:) =  ESMF_regrid_with_weights(vtmp,nemsrc+"tripole.mx025."+vecgrids(1)+".to.Ct.bilinear.nc",False)
      end do
     end if

     delete(uvar)
     delete(vvar)
     delete(ut)
     delete(vt)
     ;delete(uvmask)
    ; nv
    ;end do
   ; ii
   ;end do
  ; jj
  ;end do
end if
exit
end
