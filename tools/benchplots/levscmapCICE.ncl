undef("levscmapCICE")

function levscmapCICE(name:string,hemisphere:string)
begin

         levs = new((/3/),float)
         cmap = new((/1/),string)
    latent_levs=(/ -12.00,    12.00,  0.500/)
     ; sens_levs=(/ -20.00,    20.00,  1.000/) 
      sens_levs=(/   0.00,    40.00,  1.000/) 
      ;evap_levs=(/  -0.2,     0.2,  0.02/) 
      evap_levs=(/  -0.5,     0.5,  0.05/) 
      lwup_levs=(/-300.00,     0.00, 10.000/)
      dswf_levs=(/ 230.00,   380.00,  5.000/)
      dlwf_levs=(/ 200.00,   340.00,  5.000/)
   topmelt_levs=(/   0.00,     5.00,  0.200/)
   ;botmelt_levs=(/   0.00,     3.00,  0.100/)
   botmelt_levs=(/   0.00,     3.00,  0.050/)
      tsfc_levs=(/ -30.00,     0.00,  1.000/)
   ;    sst_levs=(/   0.00,    16.00,  0.500/)
       sst_levs=(/  -1.80,     0.00,  0.050/)
       sss_levs=(/ 28.40,  35.20,  0.200/)
      snow_levs=(/   0.00,     0.50,  0.010/)
      tair_levs=(/ -20.00,    10.00,  0.500/)
       alb_levs=(/  20.00,    80.00,  2.000/)
       ice_levs=(/   0.00,     1.00,  0.100/)
   ;    ice_levs=(/   0.90,     1.00,  0.005/)
   ;     hi_levs=(/   0.00,     4.00,  0.500/)
   if(hemisphere .eq. "SH")then
        hi_levs=(/   0.00,     2.00,  0.100/)
   else
        hi_levs=(/   0.00,     4.00,  0.200/)
   end if
     fpond_levs=(/   0.00,     1.00,  0.100/)
     hpond_levs=(/   0.00,     0.50,  0.050/)
   ;  fpond_levs=(/   0.00,     0.50,  0.050/)
   ;  hpond_levs=(/   0.00,     0.20,  0.020/)
     frzmlt_levs=(/ -20.00,  20.00,  2.000/)
     freeze_levs=(/   0.00,   2.00,  0.200/)
       tair_levs=(/  -2.00,  10.00,  0.500/)

   ; default
    cmap = "MPL_rainbow"

   if(name .eq. "evap_ai")then
    levs = evap_levs
    cmap = "BlWhRe"
   end if
   if(name .eq. "flat_ai")then
    levs = latent_levs
    cmap = "BlWhRe"
   end if
   if(name .eq. "fsens_ai")then
    levs = sens_levs
    ;cmap = "BlWhRe"
   end if
   if(name .eq. "flwup_ai")then
    levs = lwup_levs
    cmap = "MPL_rainbow"
   end if
   if(name .eq. "fswdn")then
    levs = dswf_levs
   end if
   if(name .eq. "flwdn")then
    levs = dlwf_levs
   end if
   if(name .eq. "meltt")then 
    levs = topmelt_levs
    cmap = "WhBlGrYeRe"
   end if
   if(name .eq. "meltb")then 
    levs = botmelt_levs
   end if
   if(name .eq. "Tsfc")then
    levs = tsfc_levs
   end if
   if(name .eq. "hs")then
    levs = snow_levs
    cmap = "WhBlGrYeRe"
   end if
   if(name .eq. "Tair")then
    levs = tair_levs
   end if
   if(name .eq. "albsni")then
    levs = alb_levs
    cmap = "WhBlGrYeRe"
   end if
   if(name .eq. "aice")then
    cmap = "WhBlGrYeRe"
    levs = ice_levs
   end if
   if(name .eq. "hi")then
    cmap = "WhBlGrYeRe"
    levs = hi_levs
   end if
   if(name .eq. "apond_ai" .or. name .eq. "apond")then
    cmap = "WhBlGrYeRe"
    levs = fpond_levs
   end if
   if(name .eq. "hpond_ai" .or. name .eq. "hpond")then
    cmap = "WhBlGrYeRe"
    levs = hpond_levs
   end if
   if(name .eq. "sst")then
    levs = sst_levs
   end if
   if(name .eq. "sss")then
    levs = sss_levs
   end if
   if(name .eq. "frzmlt")then
    levs = frzmlt_levs
    cmap = "ncl_default"
   end if
   if(name .eq. "congel" .or. name .eq. "frazil")then
    levs = freeze_levs
    cmap = "ncl_default"
   end if
   if(name .eq. "tair")then
    levs = tair_levs
    cmap = "ncl_default"
   end if

   ;print(name+"  "+levs)
   ;print(name+"  "+cmap)
   return( [/levs,cmap/])
end
