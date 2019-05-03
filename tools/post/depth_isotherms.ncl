begin

   ; ocean model output location 
        dirsrc = "/scratch3/NCEPDEV/stmp2/Denise.Worthen/BM1_ocn/"

          date = "20120101"

   ocnfilelist = systemfunc("ls "+dirsrc+"gfs."+idate+"/00/"+"ocn*.nc")
          ocnf = addfiles(ocnfilelist,"r")

      temp = ocnf[0]->temp
       z_l = ocnf[0]->z_l

   cosrot = ocnf[0]->cosrot
   sinrot = ocnf[0]->sinrot



  isoTlevs = (/2.5, 5.0, 10.0, 15.0 20.0, 25.0, 28.0/)
    nisoTd = dimsizes(isoTlevs)
  ;isoTlevs = ocnvar@units
  ;isoTlevs!0 = "isoTlevs"

 ; depth = conform(ocnvar, z_l, 1)

end
  copy_VarCoords(ocnvar, depth)
  printVarSummary(depth)
  isot = int2p_n_Wrap(ocnvar,depth,isoTlevs,0,1)
  printVarSummary(isot)

