;************************************************
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 

  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_inv_string.ncl"

  load "levscmapCICE.ncl"
;************************************************
 
  if(hemi .eq. "NH")then
    minlat = 66.0
  else
    maxlat = -50.0
  end if

  ; these are in degrees
     rfactor=1.0
   tlatname="TLAT"
   tlonname="TLON"
   ulatname="ULAT"
   ulonname="ULON"

  ; make faster for debugging
   iskp = 1
   jskp = 1

begin
;************************************************
; make a list of the files from the run
; this will 'catch' the initial condition file
; until Bin rewrites the workflow to give the 
; initial condition file a more unique name
;************************************************

   ;using ice2 will prevent the ic file from being included
   icefilelist = systemfunc("ls "+dirname+"ice2*.nc")
          icef = addfiles(icefilelist,"r")
        ntimes = dimsizes(icefilelist)
    ;print("ntimes = "+ntimes)
    ;print(icefilelist)

    ; get a string index
          i1 = str_index_of_substr(icefilelist(0), ".nc", 0)
       idate = icefilelist(0)
    initdate = str_get_cols(idate, i1-11, i1-1)
   ;print(initdate)

   ; use the first file to retrieve grid variables
    tmask = icef[0]->tmask(::jskp,::iskp)
    tarea = icef[0]->tarea(::jskp,::iskp)
    uarea = icef[0]->uarea(::jskp,::iskp)
  ; in radians, on t-grid
   rotang = icef[0]->ANGLET(::jskp,::iskp)

         times = icef[:]->time
     timeUnits = icef[0]->time@units
  timeCalendar = icef[0]->time@calendar

      nvars = dimsizes(varnames)
    runname = "benchmark"
   wkstitle = runname

;************************************************
; add a poly line showing the polar seem along
; great circle path from (lon1,lat1) to (lon2,lat2)
;************************************************

   lon1=  60.00
   lon2=-120.00
   lat1=  60.00
   lat2=  60.00

;************************************************
; define the map
;************************************************

  wkstype="png"
  wkstype@wkWidth=800
  wkstype@wkHeight=800
  wks   = gsn_open_wks (wkstype,wkstitle)

  npanels = dimsizes(varnames)
  ; maybe not all used
    ncols = 2
    nrows = 2

     plot = new(npanels,graphic)
      dum = new(npanels,graphic)

  res                      = True        ; plot mods desired
  res@gsnDraw              = False       ; do not draw until the end
  res@gsnFrame             = False       ; do not automatically advance 'frame'

  res@cnFillMode = "MeshFill"
  ;res@cnFillMode = "RasterFill"
  res@mpFillOn   = True

  res@cnFillOn           = True                  ; color fill
  res@cnLinesOn          = False                 ; no contour lines
  res@cnLineLabelsOn     = False                 ; no contour lines

  res@gsnPolar   =  hemi                         ; specify the hemisphere
  if(hemi .eq. "NH")then
   res@mpMinLatF    =  minlat                    ; minimum lat to plot
   res@mpCenterLonF = -45
  end if
  if(hemi .eq. "SH")then
   res@mpMaxLatF    = maxlat  
   res@mpCenterLonF = 0.0
  end if

  res@gsnRightString        = ""
  res@gsnCenterString       = ""

  res@gsnPolarLabelFontHeightF = 0.020           ; change font height of labels
  res@gsnStringFontHeightF     = 0.025
  res@lbLabelFontHeightF       = 0.020     

;************************************************
; Resources for the polyline
;************************************************

  res@mpGreatCircleLinesOn = True
  pres                  = True           ; polyline mods desired
  pres@gsLineThicknessF = 2.0            ; line thickness
  pres@gsLineColor      = "black"        ; color of lines
  pres@gsLineDashPattern =  1            ; dashes

;************************************************
;
;************************************************

  resP                     = True                ; modify the panel plot
   
  resP@gsnFrame            = False
  resP@gsnPanelLabelBar    = False               ; add common colorbar

  resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelXWhiteSpacePercent = 5

;************************************************
; loop over files
;************************************************

  do nt = 0,ntimes-1,4
  ;do nt = 0,12-1,4
       i1 = nt
       i2 = i1+3
    do nv = 0,nvars-1
         varname = varnames(nv)+"_h"
    ;average 4 files to make daily mean
          icevar = dim_avg_n_Wrap(icef[:]->$varname$(i1:i2,::jskp,::iskp),0)
    icevar@lon2d = icef[0]->$tlonname$(::jskp,::iskp)*rfactor
    icevar@lat2d = icef[0]->$tlatname$(::jskp,::iskp)*rfactor
     labelstring = icevar@long_name+" ("+icevar@units+")  "+"["+varname+"]"
    icevar@units = ""
    icevar@long_name = ""
    ;print("nv = "+nv+"  "+labelstring)
    ;printVarSummary(icevar)

   ; array contour levels for variables
       varname = varnames(nv)
       pltlevs = new((/3/),float)
       freturn = levscmapCICE(varname,hemi)
    pltlevs(:) = freturn[0]
    colormap   = freturn[1]

    if(varname .eq. "aice_h" .or. varname .eq. "aice")then
     res@cnLevelSelectionMode        = "ExplicitLevels"
     res@cnLevels                    = (/     0.15,0.50,0.7,0.8,0.85,0.9,0.950/)
     res@lbLabelStrings              = (/0.00,0.15,0.50,0.7,0.8,0.85,0.9,0.950/)+""
    else
     res@cnLevelSelectionMode = "ManualLevels"   ; set the contour levels with the following 3 resources
     res@cnMinLevelValF       = pltlevs(0)       ; set the minimum contour level
     res@cnMaxLevelValF       = pltlevs(1)       ; set the maximum contour level
     res@cnLevelSpacingF      = pltlevs(2)       ; set the interval between contours
    end if

     res@gsnLeftString  = labelstring
     res@gsnLeftStringFontHeightF = 0.025

    plot(nv) = gsn_csm_contour_map_polar(wks,icevar,res)    ; create the plot
     dum(nv) = gsn_add_polyline(wks,plot(nv),(/lon1, lon2/) , (/lat1, lat2/) ,pres)

     if(varname .eq. "aice_h" .or. varname .eq. "aice")
      delete(res@cnLevels)
      delete(res@lbLabelStrings)
     end if
    end do
    resP@gsnPanelMainString = cd_string(times(nt),"%Y%N%D")
    gsn_panel(wks,plot,(/nrows,ncols/),resP)
    frame(wks)
  end do

;************************************************
; make a movie
;************************************************

   imcmd = "convert -delay 25 "+wkstitle+".*png "+outname+initdate+".gif"
   system(imcmd)
   system("rm *.png")

exit
end
