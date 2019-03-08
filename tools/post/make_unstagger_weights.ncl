load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
    ;grdsrc  = "/scratch3/NCEPDEV/marine/noscrub/Denise.Worthen/EMC_CICE/tools/grid/"
    ;grdfile = "grid_cice_NEMS_mx025.nc"
    grdsrc  = "/scratch3/NCEPDEV/stmp2/Denise.Worthen/BenchPostTest/"
    grdfile = "ice2013070106.01.2013070100.nc"
    print(grdsrc+grdfile)

         method = "conserve"
    srcFileName = grdsrc+grdfile
    dstFileName = grdsrc+grdfile
   
    gf = addfile(grdsrc+grdfile,"r")
     lonu = gf->ULON
     latu = gf->ULAT
    uarea = gf->uarea
     lont = gf->TLON
     latt = gf->TLAT
    tarea = gf->tarea

    wgtFile = "Bu2Ct."+method+".nc"

    Opt                      = True
    Opt@InterpMethod         = method
    ;Opt@ForceOverwrite       = True
    ;Opt@PrintTimings         = True
    Opt@SrcGridArea  = uarea
    Opt@DstGridArea  = tarea
           
    Opt@SrcGridLat = latu
    Opt@SrcGridLon = lonu

    Opt@DstGridLat = latt
    Opt@DstGridLon = lont

    ESMF_regrid_gen_weights(srcFileName,dstFileName,wgtFile,Opt)
end
