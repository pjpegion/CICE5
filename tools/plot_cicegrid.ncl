;************************************************
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 

  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
;************************************************
begin
;************************************************
; read in netCDF file
;************************************************

  ; rad 2 degree
  rad2deg = 180.0d0/acos(-1.0d0)

  ;grres = "mx1"
  grres = "mx025"

  varname = "htn"
  ;varname = "hte"

;************************************************
; read in netCDF file
;************************************************

     rtsrc = "$WORK/EMC_CICE/tools/"
  ; get ang from cice grid file 
   varsrc  = "grid_cice_NEMS_"+grres+".nc"

  ; must be run in debug mode so that Cu,Cv grids
  ; are available
  ; lat,lon are degrees
  ; htn, hte are stored in cm but will be converted to km
  ; for plotting
      htnname = "htn"
   htnlonname = "lonCv"
   htnlatname = "latCv"
      htnlevs = (/6.4, 6.6, 0.005/) 

     htename  = "hte"
   htelonname = "lonCu"
   htelatname = "latCu"
      htelevs = (/2.0, 26.0, 1.00/) 

     filesrc = addfile(rtsrc+varsrc,"r")
   plottitle = rtsrc+varsrc

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
   lat1 = 30.00
   lat2 = 30.00

;************************************************
; define the map
;************************************************
  wkstype="X11"
  ;wkstype="png"
  wkstype@wkWidth=800
  wkstype@wkHeight=800
  wks   = gsn_open_wks (wkstype,"grid")              ; send graphics to X11 window

  npanels = 1
     plot = new(npanels,graphic)
      dum = new(npanels,graphic)

  res                      = True        ; plot mods desired
  res@gsnDraw              = False       ; do not draw until the end
  res@gsnFrame             = False       ; do not automatically advance 'frame'

  ;res@cnFillMode = "RasterFill"
  res@mpFillOn   = True

  if(varname .eq. "htn")then
   res@gsnPolar   = "NH"                          ; specify the hemisphere
   res@mpMinLatF  =  88                           ; minimum lat to plot
   res@mpCenterLonF = -45
   levs = htnlevs
  end if
  if(varname .eq. "hte")then
   res@mpMinLatF            =  30.0      ; range to zoom in on
   res@mpMaxLatF            =  78.0
   res@mpMinLonF            =  50.0
   res@mpMaxLonF            =  70.0
   levs = htelevs
  end if

  res@cnFillDrawOrder   = "PreDraw"             ; draw contours first
  res@cnFillOn          = False                  ; color fill
  res@cnLinesOn         = True                 ; no contour lines
  res@cnLineLabelsOn    = False

  res@cnLevelSelectionMode = "ManualLevels"     ; set the contour levels with the following 3 resources
  res@cnMinLevelValF  = levs(0)                      ; set the minimum contour level
  res@cnMaxLevelValF  = levs(1)                      ; set the maximum contour level
  res@cnLevelSpacingF = levs(2)                      ; set the interval between contours

  res@cnFillPalette        = "NCV_jet"              ; choose colormap

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

         ii = 0
  if(varname .eq. "htn")then
         var = filesrc->$htnname$(:,:)
      ; convert to km
         var = var/1.0d5
   var@lon2d = filesrc->$htnlonname$
   var@lat2d = filesrc->$htnlatname$
  var@units     = " "

  plot(ii) = gsn_csm_contour_map_polar(wks,var,res)    ; create the plot
   dum(ii) = gsn_add_polyline(wks,plot(ii),(/lon1, lon2/) , (/lat1, lat2/) ,pres)
  end if

  if(varname .eq. "hte")then
         var = filesrc->$htename$(:,:)
      ; convert to km
         var = var/1.0d5
   var@lon2d = filesrc->$htelonname$
   var@lat2d = filesrc->$htelatname$
  var@units     = " "

  plot(ii) = gsn_csm_contour_map(wks,var,res)    ; create the plot
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
