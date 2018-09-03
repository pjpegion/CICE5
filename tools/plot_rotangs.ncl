;************************************************
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 

  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
;************************************************
begin
;************************************************
;
;************************************************

  ; rad 2 degree
  rad2deg = 180.0d0/acos(-1.0d0)
 
  grres = "mx1"
  ;grres = "mx025"

  ;angle2plot = "angle"
  ; only after debug mode
  angle2plot = "anglet"

;************************************************
; read in netCDF file
;************************************************

     rtsrc = "$WORK/EMC_CICE/tools/"
  ; get ang from cice grid file 
   varsrc  = "grid_cice_NEMS_"+grres+".nc"

      angname = "angle"
   anglonname = "ulon"
   anglatname = "ulat"

  ; lat,lon are degrees
     angtname  = "anglet"
   angtlonname = "lonT"
   angtlatname = "latT"

     filesrc = addfile(rtsrc+varsrc,"r")
   plottitle = rtsrc+varsrc
       levs = (/-1.6, 1.6, 0.1/)
;       levs = (/-0.6, 0.6, 0.02/)

;************************************************
; add a poly line showing the polar seem along
; great circle path from (lon1,lat1) to (lon2,lat2)
;************************************************

  if(grres .eq. "mx1")then
   lon1=  80.00
   lon2=-100.00
  else
   lon1=  60.00
   lon2=-120.00
  end if
   lat1 = 60.00
   lat2 = 60.00

;************************************************
; define the map
;************************************************

  wkstype="X11"
  ;wkstype="png"
  wkstype@wkWidth=800
  wkstype@wkHeight=800
  wks   = gsn_open_wks (wkstype,"polar")              ; send graphics to X11 window

  npanels = 1
     plot = new(npanels,graphic)
      dum = new(npanels,graphic)

  res                      = True        ; plot mods desired
  res@gsnDraw              = False       ; do not draw until the end
  res@gsnFrame             = False       ; do not automatically advance 'frame'
  res@mpFillOn             = False
  ;res@cnFillMode           = "MeshFill"
  res@cnFillMode           = "RasterFill"

  res@cnFillDrawOrder    = "PreDraw"             ; draw contours first
  res@cnFillOn           = True              ; color fill
  res@cnLinesOn          = False             ; no contour lines
  res@lbLabelBarOn       = True                 ; turn off individual colorbars
  res@lbLabelFontHeightF = 0.015

  res@gsnPolar     = "NH"                   ; specify the hemisphere
  ;res@mpMinLatF    =  66                    ; minimum lat to plot
  res@mpMinLatF    =  88                    ; minimum lat to plot
  res@mpCenterLonF = -45

  ;res@cnFillPalette   =  "BlueWhiteOrangeRed"                      ; choose colormap
  res@cnFillPalette   =  "WhBlReWh"                      ; choose colormap

  res@cnLevelSelectionMode = "ManualLevels"     ; set the contour levels with the following 3 resources
  res@cnMinLevelValF  = levs(0)                      ; set the minimum contour level
  res@cnMaxLevelValF  = levs(1)                      ; set the maximum contour level
  res@cnLevelSpacingF = levs(2)                      ; set the interval between contours

  res@gsnLeftStringFontHeightF = 0.015

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
;
;************************************************

  if(angle2plot .eq.  "angle")then
          ii = 0
          ang = filesrc->$angname$(:,:)
   ang@lon2d = filesrc->$anglonname$*rad2deg
   ang@lat2d = filesrc->$anglatname$*rad2deg
  ang@units     = " "

  plot(ii) = gsn_csm_contour_map_polar(wks,ang,res)    ; create the plot
   dum(ii) = gsn_add_polyline(wks,plot(ii),(/lon1, lon2/) , (/lat1, lat2/) ,pres)
  end if

  if(angle2plot .eq. "anglet")then
          ii = 0
         ang = filesrc->$angtname$(:,:)
   ang@lon2d = filesrc->$angtlonname$
   ang@lat2d = filesrc->$angtlatname$
  ang@units     = " "

  plot(ii) = gsn_csm_contour_map_polar(wks,ang,res)    ; create the plot
   dum(ii) = gsn_add_polyline(wks,plot(ii),(/lon1, lon2/) , (/lat1, lat2/) ,pres)
  end if

;************************************************
;
;************************************************

  resP                     = True                ; modify the panel plot

  resP@gsnPanelMainString  = plottitle
  resP@gsnPanelMainFontHeightF = 0.010            ; new resource added in NCL V6.4.0

  gsn_panel(wks,plot,(/1,1/),resP)               ; now draw as one plot

exit
end
