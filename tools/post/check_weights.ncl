; Check weights
;
; Create fields on stagger grid and see that they are 
; interpolated to Ct grid bilinearly
;
; This file still has to be loaded manually
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
;------------------------------------------------------------------
begin


;---Specify a location to use
        nemsrc     = "/scratch4/NCEPDEV/ocean/save/Denise.Worthen/NEMS_INPUT0.1/ocnicepost/"
            srcloc = nemsrc
          gridfile = "tripole.grid.mx025.nc"
       staggerlocs = (/"Ct", "Cu", "Cv", "Bu"/)

      gf = addfile(nemsrc+gridfile,"r")
   ;ii =1
   do ii = 1,dimsizes(staggerlocs)-1
         wgtfile = nemsrc+"tripole.grid.mx025."+staggerlocs(ii)+".to.Ct.bilinear.nc"
       ; get a field
         varname = "lon"+staggerlocs(ii)
            avar = gf->$varname$
       print(wgtfile+"  "+varname)
 
       avar = 0.0
       dims = dimsizes(avar)
       jmax = dims(0)
       imax = dims(1)
     
      ; will create field with integer values
      ; bilinear interpolation should result in 0.5 values 
       if(staggerlocs(ii) .eq. "Cu")then
        do i = 0,imax-1
          avar(:,i) = todouble(i)
        end do
       end if
       if(staggerlocs(ii) .eq. "Cv")then
        do j = 0,jmax-1
          avar(j,:) = todouble(j)
        end do
       end if
       if(staggerlocs(ii) .eq. "Bu")then
        do j = 0,jmax-1
         do i = 0,imax-1
          avar(j,i) = todouble(2.0*mod(i,2)*mod(j,2))
         end do
        end do
       end if

;----------------------------------------------------------------------
; Apply the weights to regrid to Ct
;----------------------------------------------------------------------

    Opt                = True
    Opt@PrintTimings   = True
    Opt@Debug          = True

    avarCt = ESMF_regrid_with_weights(avar,wgtfile,Opt)
   
    ; write to a netcdf file for examination
    outfile = nemsrc+staggerlocs(ii)+".regrid.Ct.check.nc"
    system("/bin/rm -f "+outfile)           ; remove any pre-existing file
      ncdf = addfile(outfile ,"c")          ; open output netCDF file

      dimSizes = (/jmax,imax/)
      dimNames = (/"nj", "ni"/)
      dimUnlim = (/False, False/)
      filedimdef(ncdf,dimNames,dimSizes,dimUnlim)

      filevardef(ncdf,  "avarCt", typeof(avar),(/"nj","ni"/))
      filevardef(ncdf,  "avar",   typeof(avar),(/"nj","ni"/))
      ncdf->avar      = (/avar/)
      ncdf->avarCt    = (/avarCt/)
      delete(avar)
      delete(avarCt)
   end do
exit
end

