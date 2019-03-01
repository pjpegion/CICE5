;************************************************
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 

  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;************************************************
begin

  batch = "no"

  hemi="NH"
  ;hemi="SH"

  ;minlat = 66.0
  ;minlat = 76.0
  minlat = 88.0
  ;minlat = 80.0

   rhoi = 917.0
   rhos = 330.0
;************************************************
; read in netCDF file
;************************************************

  RT1 = "/scratch3/NCEPDEV/stmp1/Denise.Worthen/"
  RT2 = "/scratch3/NCEPDEV/stmp2/Denise.Worthen/"
  RT3 = "/scratch3/NCEPDEV/stmp3/Denise.Worthen/"
  RT4 = "/scratch4/NCEPDEV/stmp4/Denise.Worthen/"

  ncat = 5

;************************************************
; read in netCDF file
;************************************************

  ;cdate = "2011040100"
  cdate = "2011090100"

   sis2 = "ice_model.res_"+cdate+".nc"

  ; this copy contains the center lats and lons in degrees, which
  ; is more convenient for plotting
  dirsrc = "/scratch4/NCEPDEV/ocean/save/Denise.Worthen/NEMS_INPUT0.1/cice_data/mx024/"
  ; load the grid file
      gr = addfile(dirsrc+"grid_cice_NEMS_mx025.nc","r")

  icelonname = "lonT"
  icelatname = "latT"

      sis2file = addfile(RT1+sis2, "r")
        icecat = sis2file->part_size(0,1:ncat,:,:)
        icethk = sis2file->h_ice(0,:,:,:)
  icecat@lon2d = gr->$icelonname$
  icecat@lat2d = gr->$icelatname$
  icethk@lon2d = gr->$icelonname$
  icethk@lat2d = gr->$icelatname$

  icethk@_FillValue = default_fillvalue("double")
  do k = 0,ncat-1
   icethk(k,:,:) = where(icethk@lat2d .ge. minlat, icecat(k,:,:)*icethk(k,:,:)/rhoi, icethk@_FillValue)
  end do

  wkstitle="restarts"+cdate
;************************************************
;
;************************************************

  names = (/"part_size","h_ice"/)

  ;ailevs = (/0.98, 1.0, 0.001/)
  ailevs = (/0.89, 0.93, 0.001/)
  ;ailevs = (/0.00, 1.0, 0.05/)
  ;hilevs = (/0.00, 4.0, 0.10/)
  ;hilevs = (/0.00, 3.0, 0.05/)
  ;hilevs = (/1.50, 1.68, 0.005/)
  hilevs = (/0.55, 0.90, 0.005/)

;************************************************
; add a poly line showing the polar seem along
; great circle path from (lon1,lat1) to (lon2,lat2)
;************************************************

  lon1=  60.00
  lon2=-120.00
  lat1 = 60.00
  lat2 = 60.00

;************************************************
; define the map
;************************************************

  if(batch .eq. "no")then
   wkstype="X11"
  else
   wkstype="png"
  end if
   ;wkstype="png"
  wkstype@wkWidth=800
  wkstype@wkHeight=800
  wks   = gsn_open_wks (wkstype,wkstitle)

  ; panels
  npanels = 2
  plot = new(npanels,graphic)
   dum = new(npanels,graphic)

  res                      = True        ; plot mods desired
  res@gsnDraw              = False       ; do not draw until the end
  res@gsnFrame             = False       ; do not automatically advance 'frame'

  res@mpFillOn   = True
  ;res@cnFillMode           = "MeshFill"
  res@cnFillMode           = "RasterFill"
  res@cnFillDrawOrder      = "PreDraw"             ; draw contours first

  res@cnFillOn           = True                  ; color fill
  res@cnLinesOn          = False                 ; no contour lines
  res@lbLabelBarOn       = True                 ; turn off individual colorbars
  res@lbLabelFontHeightF = 0.015

  res@gsnStringFontHeightF = 0.05

  res@gsnPolar   = "NH"                          ; specify the hemisphere
  res@mpMinLatF    =  minlat                          ; minimum lat to plot
  res@mpCenterLonF = -45

  res@cnFillPalette        = "MPL_rainbow"      ; choose colormap
  
  res@gsnStringFontHeightF = 0.020
;************************************************
; Resources for the polyline
;************************************************

  res@mpGreatCircleLinesOn = True
  pres                  = True           ; polyline mods desired
  pres@gsLineThicknessF = 2.0            ; line thickness
  pres@gsLineColor      = "black"         ; color of lines
  pres@gsLineDashPattern =  1            ; dashes
  ;pres@tfPolyDrawOrder = "postDraw"

;************************************************
; the subplot call
;************************************************
function subplot(ii,fname,name,pltvar,pltlevs)
begin

   res@cnLevelSelectionMode = "ManualLevels"        ; set the contour levels with the following 3 resources
   res@cnMinLevelValF       = pltlevs(0)            ; set the minimum contour level
   res@cnMaxLevelValF       = pltlevs(1)            ; set the maximum contour level
   res@cnLevelSpacingF      = pltlevs(2)            ; set the interval between contours

   ;res@gsnLeftString = name
   res@gsnCenterString = fname

   splot = gsn_csm_contour_map_polar(wks,pltvar,res)    ; create the plot

   return(splot)
end
;************************************************
;
;************************************************

  print("plotting variable "+names)

      ii = 0
  ; k=0 is open water, cats 1:5 are ice
    levs = ailevs
    name = names(ii)
  ; to plot the sum over the categories
  plotvar = dim_sum_n_Wrap(icecat,0)
  plotvar@units = " "
  plotvar@lon2d = gr->$icelonname$
  plot(ii) = subplot(ii,sis2,name,plotvar,levs)
   dum(ii) = gsn_add_polyline(wks,plot(ii),(/lon1, lon2/) , (/lat1, lat2/) ,pres)

      ii = 1
    levs = hilevs
    name = names(ii)
  ; to plot the sum over the categories
  plotvar = dim_sum_n_Wrap(icethk,0)
  plotvar@long_name = " "
  plotvar@units = " "
  plotvar@lon2d = gr->$icelonname$
  plotvar@lat2d = gr->$icelatname$
  plot(ii) = subplot(ii,sis2,name,plotvar,levs)
   dum(ii) = gsn_add_polyline(wks,plot(ii),(/lon1, lon2/) , (/lat1, lat2/) ,pres)

;************************************************
; panel the plots
;************************************************

  resP                     = True                ; modify the panel plot
                                                 ; new resource added in NCL V6.4.0
  ;resP@gsnPanelLabelBar    = True                ; add common colorbar
  ;resP@lbLabelFontHeightF  = 0.007               ; make labels smaller

 ;resP@gsnPanelYWhiteSpacePercent = 5
 ;resP@gsnPanelXWhiteSpacePercent = 5

  gsn_panel(wks,plot,(/1,2/),resP)               ; now draw as one plot

exit
end

