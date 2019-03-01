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

  minlat = 66.0
  ;minlat = 76.0
  ;minlat = 88.0
  ;minlat = 80.0

  ; define something approximating cryosat's 
  ; awful colormap

  cryosat2 = (/ (/ 45,  51,  34/), (/ 30,  47, 118/), \
                (/  6,  85, 121/), (/  0, 122,  88/), \
                (/  0, 160,  55/), (/ 78, 137,  39/), \
                (/198, 201,  42/), (/247, 219, 196/) /)/255.

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

  cice5 = "cice5_model_0.25.res_"+cdate+".nc"
  cfsv2 = "cfsv2_cice5_model_0.25.res_"+cdate+".nc"

  ; this copy contains the center lats and lons in degrees, which
  ; is more convenient for plotting
  dirsrc = "/scratch4/NCEPDEV/ocean/save/Denise.Worthen/NEMS_INPUT0.1/cice_data/mx024/"
  ; load the grid file
      gr = addfile(dirsrc+"grid_cice_NEMS_mx025.nc","r")

  icelonname = "lonT"
  icelatname = "latT"

  ; to plot ice thickness
  ; varname = "vicen"
  ; contour levels from 0,4 by 0.2
  ;   levs = (/0.0, 4.0, 0.2/)

  cice5file = addfile(RT1+cice5,"r")
  cfsv2file = addfile(RT1+cfsv2,"r")
  ;print(sis2)
  ;print(cice5)
  ;print(cfsv2)

  icefilelist = (/RT1+cice5, RT1+cfsv2/)
     icefiles = addfiles(icefilelist,"r")
       nfiles = dimsizes(icefilelist)
 print(icefilelist)

  wkstitle="restarts"+cdate
;************************************************
;
;************************************************

  rhoi = 917.0
  rhos = 330.0
  ainames = (/"aicen", "aicen"/)
  hinames = (/"vicen", "vicen"/)

  ;ailevs = (/0.85, 1.0, 0.01/)
  ailevs = (/0.00, 1.0, 0.05/)
  ;hilevs = (/0.55, 5.00, 0.0275/)
  hilevs = (/0.00, 5.0, 0.0275/)

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
  ; wkstype="png"
  wkstype@wkWidth=800
  wkstype@wkHeight=800
  wks   = gsn_open_wks (wkstype,wkstitle)

  ; panels
  npanels = nfiles
  plot = new(npanels,graphic)
   dum = new(npanels,graphic)

  res                      = True        ; plot mods desired
  res@gsnDraw              = False       ; do not draw until the end
  res@gsnFrame             = False       ; do not automatically advance 'frame'

  res@mpFillOn   = True
  res@cnFillMode           = "MeshFill"
  ;res@cnFillMode           = "RasterFill"
  res@cnFillDrawOrder      = "PreDraw"             ; draw contours first

  res@cnFillOn           = True                  ; color fill
  res@cnLinesOn          = False                 ; no contour lines
  res@lbLabelBarOn       = True                 ; turn off individual colorbars
  res@lbLabelFontHeightF = 0.015

  res@gsnStringFontHeightF = 0.05

  res@gsnPolar   = "NH"                          ; specify the hemisphere
  res@mpMinLatF    =  minlat                          ; minimum lat to plot
  res@mpCenterLonF = -45

  res@lbBoxSeparatorLinesOn = "False"
  
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
   if(name .eq. "vicen")then
    res@cnFillPalette        = cryosat2      ; choose colormap
   else
    res@cnFillPalette        = "MPL_rainbow"      ; choose colormap
   end if
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
  ;  levs = ailevs
  ; names = ainames
    levs = hilevs
   names = hinames
  print("plotting variable "+names)

      ii = 0
    name = names(ii)
  icevar = icefiles[ii]->$name$(:,:,:)
  icevar@long_name = " "
  icevar@units = " "

  ; to plot the sum over the categories
  plotvar = dim_sum_n(icevar,0)
  plotvar = where(plotvar .eq. 0.0, default_fillvalue("double"), plotvar)
  plotvar@lon2d = gr->$icelonname$
  plotvar@lat2d = gr->$icelatname$
  plot(ii) = subplot(ii,cice5,name,plotvar,levs)
   dum(ii) = gsn_add_polyline(wks,plot(ii),(/lon1, lon2/) , (/lat1, lat2/) ,pres)
  delete(icevar)

      ii = 1
    name = names(ii)
  icevar = icefiles[ii]->$name$(:,:,:)
  icevar@long_name = " "
  icevar@units = " "

  ; to plot the sum over the categories
  plotvar = dim_sum_n(icevar,0)
  plotvar = where(plotvar .eq. 0.0, default_fillvalue("double"), plotvar)
  plotvar@lon2d = gr->$icelonname$
  plotvar@lat2d = gr->$icelatname$
  plot(ii) = subplot(ii,cfsv2,name,plotvar,levs)
   dum(ii) = gsn_add_polyline(wks,plot(ii),(/lon1, lon2/) , (/lat1, lat2/) ,pres)
  delete(icevar)

;************************************************
; panel the plots
;************************************************

  resP                     = True                ; modify the panel plot
                                                 ; new resource added in NCL V6.4.0
  ;resP@gsnPanelLabelBar    = True                ; add common colorbar
  ;resP@lbLabelFontHeightF  = 0.007               ; make labels smaller

 ;resP@gsnPanelYWhiteSpacePercent = 5
 ;resP@gsnPanelXWhiteSpacePercent = 5

  gsn_panel(wks,plot,(/2,1/),resP)               ; now draw as one plot

exit
end

