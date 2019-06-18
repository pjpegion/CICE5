load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;begin

testdir="/scratch4/NCEPDEV/stmp3/Denise.Worthen/CPC/"
testdates=(/"20160201","20160215","20160601","20120801","20120801"/)
ndates=dimsizes(testdates)
print(ndates)
print("Entered convert_cice5_025.ncl")

ftmask = addfile("kmtu_cice_NEMS_mx025.nc","r")
kmt = ftmask->kmt
kmt5 = new((/5,1080,1440/),"double")
do k = 0,4
  kmt5(k,:,:) = kmt(:,:)
end do

; loop over dates
;do ii = 0,ndates-1
  cdate = testdates(0)+"00"
  fin = addfile(testdir+"cpc."+cdate+"/ice_model.res.nc","r")
; fin  = addfile("ice_model.res.nc","r")

;print(cdate)
;finname="ice_model.res_"+cdate+".nc"
;print(finname)
;fin  = addfile(finname,"r")


Opt = True
Opt@PrintTimings = True

part_size1 = fin->part_size
ndims1 = dimsizes(part_size1)
nilyr1 = 2
nilyr2 = 7
nslyr = 1
ncat = ndims1(1) - 1
nj1 = ndims1(2)
ni1 = ndims1(3)
ni2 = 1440
nj2 = 1080

h_ice1 = fin->h_ice
h_sno1 = fin->h_snow
t_surf1 = fin->t_surf
t_ice1 = fin->t_ice1
t_ice2 = fin->t_ice2

part_size = new((/ncat,nj2,ni2/),"double")
h_ice = new((/ncat,nj2,ni2/),"double")
h_sno = new((/ncat,nj2,ni2/),"double")
t_surf = new((/ncat,nj2,ni2/),"double")
t1 = new((/ncat,nj2,ni2/),"double")
t2 = new((/ncat,nj2,ni2/),"double")

do n=0,ncat-1
   part_size(n,:,:) = ESMF_regrid_with_weights(part_size1(0,n+1,:,:),"MOM0.5_to_MOM0.25_map.nc",Opt)
   h_ice(n,:,:) = ESMF_regrid_with_weights(h_ice1(0,n,:,:),"MOM0.5_to_MOM0.25_map.nc",Opt)
   h_sno(n,:,:) = ESMF_regrid_with_weights(h_sno1(0,n,:,:),"MOM0.5_to_MOM0.25_map.nc",Opt)
   t_surf(n,:,:) = ESMF_regrid_with_weights(t_surf1(0,n+1,:,:),"MOM0.5_to_MOM0.25_map.nc",Opt)
   t1(n,:,:) = ESMF_regrid_with_weights(t_ice1(0,n,:,:),"MOM0.5_to_MOM0.25_map.nc",Opt)
   t2(n,:,:) = ESMF_regrid_with_weights(t_ice2(0,n,:,:),"MOM0.5_to_MOM0.25_map.nc",Opt)
end do

tmp = fin->u_ice
uvel = ESMF_regrid_with_weights(tmp,"MOM0.5_to_MOM0.25_map.nc",Opt)
tmp = fin->v_ice
vvel = ESMF_regrid_with_weights(tmp,"MOM0.5_to_MOM0.25_map.nc",Opt)
tmp = fin->sig11
sig11 = ESMF_regrid_with_weights(tmp,"MOM0.5_to_MOM0.25_map.nc",Opt)
tmp = fin->sig22
sig22 = ESMF_regrid_with_weights(tmp,"MOM0.5_to_MOM0.25_map.nc",Opt)
tmp = fin->sig12
sig12 = ESMF_regrid_with_weights(tmp,"MOM0.5_to_MOM0.25_map.nc",Opt)

; tmp changes dimensions
delete(tmp)
tmp = fin->flux_sw_vis_dir
swvdr = ESMF_regrid_with_weights(tmp,"MOM0.5_to_MOM0.25_map.nc",Opt)
tmp = fin->flux_sw_vis_dif
swvdf = ESMF_regrid_with_weights(tmp,"MOM0.5_to_MOM0.25_map.nc",Opt)
tmp = fin->flux_sw_nir_dir
swidr = ESMF_regrid_with_weights(tmp,"MOM0.5_to_MOM0.25_map.nc",Opt)
tmp = fin->flux_sw_nir_dif
swidf = ESMF_regrid_with_weights(tmp,"MOM0.5_to_MOM0.25_map.nc",Opt)
; remove negative values (present in cpc)
swvdr = where(swvdr .lt. 0.0, 0.0, swvdr)
swvdf = where(swvdf .lt. 0.0, 0.0, swvdf)
swidr = where(swidr .lt. 0.0, 0.0, swidr)
swidf = where(swidf .lt. 0.0, 0.0, swidf)
; These need to be created

; We don't have frz_onst from CFS runs
frz_onset = uvel*0.d0

aicen = new((/ncat,nj2,ni2/),double)
vicen = new((/ncat,nj2,ni2/),double)
vsnon = new((/ncat,nj2,ni2/),double)
Tsfcn = new((/ncat,nj2,ni2/),double)

aicen = where(ismissing(part_size),0.0d0,part_size)
print(aicen(:,74:75,694:695))
print(kmt5(:,74:75,694:695))
; aicen = where(ismissing(kmt5),0.0d0,aicen)
do k=0,4
do j=0,1079
do i=0,1439
 aicen(k,j,i) = where(kmt5(k,j,i) < 0.99, aicen(k,j,i), 0)
end do
end do
end do

print(aicen(:,74:75,694:695))

vicen = where(ismissing(h_ice),0.0d0,h_ice*aicen)
vsnon = where(ismissing(h_sno),0.0d0,h_sno*aicen)
Tsfcn = where(ismissing(t_surf),0.0d0,t_surf - 273.15d0)

Tsfcn = where(Tsfcn.gt.0.0d0,0.0d0,Tsfcn)


Tin = new((/nilyr2,ncat,nj2,ni2/),double)
qice = new((/nilyr2,ncat,nj2,ni2/),double)
sice = new((/nilyr2,ncat,nj2,ni2/),double)
qsno = new((/nslyr,ncat,nj2,ni2/),double)
qice = 0.d0
sice = 0.d0
qsno = 0.d0

; Linearly interpolate two layers to 7.
xm = (t1-t2)*2.0d0
b = t1 - xm*0.75d0
Tin(0,:,:,:) = xm*0.929d0 + b 
Tin(1,:,:,:) = xm*0.786d0 + b
Tin(2,:,:,:) = xm*0.643d0 + b
Tin(3,:,:,:) = xm*0.500d0 + b
Tin(4,:,:,:) = xm*0.357d0 + b
Tin(5,:,:,:) = xm*0.214d0 + b
Tin(6,:,:,:) = xm*0.071d0 + b

Tin = where(Tin.lt.0.,Tin,0.)
Tin = where(Tin.ge.-100.,Tin,0.)

; Salinity
saltmax = 3.2d0
nsal = 0.407d0
msal = 0.573d0
pi = atan(1.0d0)*4.0d0

salinz = new((/nilyr2/),double)
do k=0,nilyr2-1
   zn = (int2dble(k+1)-0.5d0)/int2dble(nilyr2)
   salinz(k) = (saltmax/2.d0)*(1.d0-cos(pi*zn^(nsal/(msal+zn))))
   sice(k,:,:,:) = salinz(k)
end do

rhoi = 917.d0
rhos = 330.d0
cp_ice = 2106.d0
cp_ocn = 4218.d0
Lfresh = 3.34d5


; Convert temperature to enthalpy
Tmltz = salinz / (-18.48 + (0.01848*salinz))
print(salinz)
print(Tmltz)
do k=0,nilyr2-1
   Tin(k,:,:,:) = where(Tin(k,:,:,:).gt.Tmltz(k),Tmltz(k),Tin(k,:,:,:))
;  qice(k,:,:,:) = -rhoi*(cp_ice*(Tmltz(k)-Tin(k,:,:,:)) + Lfresh*(1.0d0-Tmltz(k)/Tin(k,:,:,:)) - cp_ocn*Tmltz(k)) 
   qice(k,:,:,:) = rhoi*cp_ice*Tin(k,:,:,:) - rhoi*Lfresh
   qice(k,:,:,:) = where(vicen.gt.0.d0,qice(k,:,:,:),0.)
end do

qice = where(ismissing(qice),0.0d0,qice)

do k=0,nslyr-1
   qsno(k,:,:,:) = -rhos*(Lfresh - cp_ice*Tsfcn(:,:,:))
   qsno(k,:,:,:) = where(vsnon.gt.0.d0,qsno(k,:,:,:),0.)
end do

print(Tsfcn(:,400,163))
print(Tin(:,:,400,163))
print(qsno(0,:,400,163))
print(qice(0,:,400,163))

stressp = sig11(0,:,:)+sig22(0,:,:)
stressm = sig11(0,:,:)-sig22(0,:,:)
stress12 = sig12(0,:,:)

stressp_1 = stressp*0.
stressp_2 = stressp*0.
stressp_3 = stressp*0.
stressp_4 = stressp*0.
stressm_1 = stressp*0.
stressm_2 = stressp*0.
stressm_3 = stressp*0.
stressm_4 = stressp*0.
stress12_1 = stressp*0.
stress12_2 = stressp*0.
stress12_3 = stressp*0.
stress12_4 = stressp*0.

x1 = flt2dble(fspan(1,ni2,ni2))
x1p = x1+0.5
x1m = x1-0.5
y1 = flt2dble(fspan(1,nj2,nj2))

y1p = y1+0.5
y1m = y1-0.5

print ("starting stressp_1 linint2")
stressp_1 = linint2(x1,y1,stressp,True,x1p,y1p,0)
print ("starting stressp_2 linint2")
stressp_2 = linint2(x1,y1,stressp,True,x1m,y1p,0)
print ("starting stressp_3 linint2")
stressp_3 = linint2(x1,y1,stressp,True,x1m,y1m,0)
print ("starting stressp_4 linint2")
stressp_4 = linint2(x1,y1,stressp,True,x1p,y1m,0)
print ("starting stressm_1 linint2")
stressm_1 = linint2(x1,y1,stressm,True,x1p,y1p,0)
print ("starting stressm_2 linint2")
stressm_2 = linint2(x1,y1,stressm,True,x1m,y1p,0)
print ("starting stressm_3 linint2")
stressm_3 = linint2(x1,y1,stressm,True,x1m,y1m,0)
print ("starting stressm_4 linint2")
stressm_4 = linint2(x1,y1,stressm,True,x1p,y1m,0)
print ("starting stress12_1 linint2")
stress12_1 = linint2(x1,y1,stress12,True,x1p,y1p,0)
print ("starting stress12_2 linint2")
stress12_2 = linint2(x1,y1,stress12,True,x1m,y1p,0)
print ("starting stress12_2 linint2")
stress12_3 = linint2(x1,y1,stress12,True,x1m,y1m,0)
print ("starting stress12_4 linint2")
stress12_4 = linint2(x1,y1,stress12,True,x1p,y1m,0)


stressp_1 = where(ismissing(stressp_1),0.0d0,stressp_1)
stressp_2 = where(ismissing(stressp_2),0.0d0,stressp_2)
stressp_3 = where(ismissing(stressp_3),0.0d0,stressp_3)
stressp_4 = where(ismissing(stressp_4),0.0d0,stressp_4)
stressm_1 = where(ismissing(stressm_1),0.0d0,stressm_1)
stressm_2 = where(ismissing(stressm_2),0.0d0,stressm_2)
stressm_3 = where(ismissing(stressm_3),0.0d0,stressm_3)
stressm_4 = where(ismissing(stressm_4),0.0d0,stressm_4)
stress12_1 = where(ismissing(stress12_1),0.0d0,stress12_1)
stress12_2 = where(ismissing(stress12_2),0.0d0,stress12_2)
stress12_3 = where(ismissing(stress12_3),0.0d0,stress12_3)
stress12_4 = where(ismissing(stress12_4),0.0d0,stress12_4)

aice = dim_sum_n(aicen,0)
iceumask = aice*0.d0
iceumask = where(aice.gt.1.0d-11,1.,0.)

; Initialize these to zero
coszen = uvel*0.
scale_factor = uvel*0.
;swvdr = uvel*0.
;swvdf = uvel*0.
;swidr = uvel*0.
;swidf = uvel*0.
strocnxT = uvel*0.
strocnyT = uvel*0.

file_atts = getvaratts(fin)
natts = dimsizes(file_atts)

setfileoption("nc","format","LargeFile")

; fout = addfile("cice5_model_0.25.res.nc","c")

foutname = testdir+"cpc."+cdate+"/cice5_model_0.25.res_"+cdate+".nc"
;foutname="cice5_model_0.25.res_"+cdate+".nc"
system("/bin/rm -f "+foutname)
fout = addfile(foutname,"c")
print(foutname)


setfileoption(fout,"DefineMode",True)

do iatt=0,natts-1
   fout@$file_atts(iatt)$ = fin@$file_atts(iatt)$
end do

dimNames = (/"nilyr","ncat","nj","ni"/)
dimSizes = (/nilyr2,ncat,nj2,ni2/)
dimUnlim = (/False,False,False,False/)
filedimdef(fout,dimNames,dimSizes,dimUnlim)

filevardef(fout,"aicen",typeof(aicen),(/"ncat","nj","ni"/))
filevardef(fout,"vicen",typeof(vicen),(/"ncat","nj","ni"/))
filevardef(fout,"vsnon",typeof(vsnon),(/"ncat","nj","ni"/))
filevardef(fout,"Tsfcn",typeof(Tsfcn),(/"ncat","nj","ni"/))
filevardef(fout,"uvel",typeof(uvel),(/"nj","ni"/))
filevardef(fout,"vvel",typeof(vvel),(/"nj","ni"/))
filevardef(fout,"scale_factor",typeof(scale_factor),(/"nj","ni"/))
filevardef(fout,"coszen",typeof(coszen),(/"nj","ni"/))
filevardef(fout,"swvdr",typeof(swvdr),(/"nj","ni"/))
filevardef(fout,"swvdf",typeof(swvdf),(/"nj","ni"/))
filevardef(fout,"swidr",typeof(swidr),(/"nj","ni"/))
filevardef(fout,"swidf",typeof(swidf),(/"nj","ni"/))
filevardef(fout,"strocnxT",typeof(strocnxT),(/"nj","ni"/))
filevardef(fout,"strocnyT",typeof(strocnyT),(/"nj","ni"/))
filevardef(fout,"stressp_1",typeof(stressp_1),(/"nj","ni"/))
filevardef(fout,"stressp_2",typeof(stressp_2),(/"nj","ni"/))
filevardef(fout,"stressp_3",typeof(stressp_3),(/"nj","ni"/))
filevardef(fout,"stressp_4",typeof(stressp_4),(/"nj","ni"/))
filevardef(fout,"stressm_1",typeof(stressm_1),(/"nj","ni"/))
filevardef(fout,"stressm_2",typeof(stressm_2),(/"nj","ni"/))
filevardef(fout,"stressm_3",typeof(stressm_3),(/"nj","ni"/))
filevardef(fout,"stressm_4",typeof(stressm_4),(/"nj","ni"/))
filevardef(fout,"stress12_1",typeof(stress12_1),(/"nj","ni"/))
filevardef(fout,"stress12_2",typeof(stress12_2),(/"nj","ni"/))
filevardef(fout,"stress12_3",typeof(stress12_3),(/"nj","ni"/))
filevardef(fout,"stress12_4",typeof(stress12_4),(/"nj","ni"/))
filevardef(fout,"iceumask",typeof(iceumask),(/"nj","ni"/))

filevardef(fout,"qice001",typeof(qice),(/"ncat","nj","ni"/))
filevardef(fout,"qice002",typeof(qice),(/"ncat","nj","ni"/))
filevardef(fout,"qice003",typeof(qice),(/"ncat","nj","ni"/))
filevardef(fout,"qice004",typeof(qice),(/"ncat","nj","ni"/))
filevardef(fout,"qice005",typeof(qice),(/"ncat","nj","ni"/))
filevardef(fout,"qice006",typeof(qice),(/"ncat","nj","ni"/))
filevardef(fout,"qice007",typeof(qice),(/"ncat","nj","ni"/))
filevardef(fout,"sice001",typeof(sice),(/"ncat","nj","ni"/))
filevardef(fout,"sice002",typeof(sice),(/"ncat","nj","ni"/))
filevardef(fout,"sice003",typeof(sice),(/"ncat","nj","ni"/))
filevardef(fout,"sice004",typeof(sice),(/"ncat","nj","ni"/))
filevardef(fout,"sice005",typeof(sice),(/"ncat","nj","ni"/))
filevardef(fout,"sice006",typeof(sice),(/"ncat","nj","ni"/))
filevardef(fout,"sice007",typeof(sice),(/"ncat","nj","ni"/))
filevardef(fout,"qsno001",typeof(qsno),(/"ncat","nj","ni"/))
filevardef(fout,"frz_onset",typeof(frz_onset),(/"nj","ni"/))
filevardef(fout,"Tin",typeof(Tin),(/"nilyr","ncat","nj","ni"/))

setfileoption(fout,"DefineMode",False)
;if(1.eq.0)then
; setfileoption(fout,"DefineMode",True)

; strtWrt = systemfunc("date")
; print (strtWrt)
print ("writing aicen")
fout->aicen = (/aicen/)
print ("writing vicen")
fout->vicen = (/vicen/)
print ("writing vsnon")
fout->vsnon = (/vsnon/)
print ("writing Tsfcn")
fout->Tsfcn = (/Tsfcn/)
print ("writing uvel")
fout->uvel = (/uvel/)
print ("writing vvel")
fout->vvel = (/vvel/)
print ("writing scale_factor")
fout->scale_factor = (/scale_factor/)
print ("writing coszen")
fout->coszen = (/coszen/)
print ("writing swvdr")
fout->swvdr = (/swvdr/)
print ("writing swvdf")
fout->swvdf = (/swvdf/)

; wallClockElapseTime(strtWrt, "write 10 records", 0)

print ("writing swidr")
fout->swidr = (/swidr/)
print ("writing swidf")
fout->swidf = (/swidf/)
print ("writing strocnxT")
fout->strocnxT = (/strocnxT/)
print ("writing strocnyT")
fout->strocnyT = (/strocnyT/)
print ("writing stressp_1")
fout->stressp_1 = (/stressp_1/)
print ("writing stressp_2")
fout->stressp_2 = (/stressp_2/)
print ("writing stressp_3")
fout->stressp_3 = (/stressp_3/)
print ("writing stressp_4")
fout->stressp_4 = (/stressp_4/)
print ("writing stressm_1")
fout->stressm_1 = (/stressm_1/)
print ("writing stressm_2")
fout->stressm_2 = (/stressm_2/)
print ("writing stressm_3")
fout->stressm_3 = (/stressm_3/)
print ("writing stressm_4")
fout->stressm_4 = (/stressm_4/)
print ("writing stress12_4")
fout->stress12_4 = (/stress12_4/)
print ("writing stress12_1")
fout->stress12_1 = (/stress12_1/)
print ("writing stress12_2")
fout->stress12_2 = (/stress12_2/)
print ("writing stress12_3")
fout->stress12_3 = (/stress12_3/)
print ("writing iceumask")
fout->iceumask = (/iceumask/)
print ("writing frz_onset")
fout->frz_onset = (/frz_onset/)
print ("writing qice001")
fout->qice001 = (/qice(0,:,:,:)/)
print ("writing qice002")
fout->qice002 = (/qice(1,:,:,:)/)
print ("writing qice003")
fout->qice003 = (/qice(2,:,:,:)/)
print ("writing qice004")
fout->qice004 = (/qice(3,:,:,:)/)
print ("writing qice005")
fout->qice005 = (/qice(4,:,:,:)/)
print ("writing qice006")
fout->qice006 = (/qice(5,:,:,:)/)
print ("writing qice007")
fout->qice007 = (/qice(6,:,:,:)/)
print ("writing Tin")
fout->Tin = (/Tin/)
print ("writing sice001")
fout->sice001 = (/sice(0,:,:,:)/)
print ("writing sice002")
fout->sice002 = (/sice(1,:,:,:)/)
print ("writing sice003")
fout->sice003 = (/sice(2,:,:,:)/)
print ("writing sice004")
fout->sice004 = (/sice(3,:,:,:)/)
print ("writing sice005")
fout->sice005 = (/sice(4,:,:,:)/)
print ("writing sice006")
fout->sice006 = (/sice(5,:,:,:)/)
print ("writing sice007")
fout->sice007 = (/sice(6,:,:,:)/)
print ("writing qsno001")
fout->qsno001 = (/qsno(0,:,:,:)/)
;end if
print ("The END")

delete([/part_size,h_ice,h_sno,t_surf,t1,t2,aicen,vicen,vsnon,Tsfcn,Tin,qice,sice,qsno,salinz/])
; date loop
;end do
;end
