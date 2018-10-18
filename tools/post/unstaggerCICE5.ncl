undef("unstagger")
; based on wrf unstagger
; CICE5 native velocities are colocated on corners, 2D only
function unstaggerCICE5(varin:numeric,areaU:numeric,areaT:numeric)
begin
  tdims = dimsizes(varin)
  ;print("varin = "+tdims(0)+"  "+tdims(1))
  varout = new((/tdims(0),tdims(1)/),float)

  dims = dimsizes(varin)
    nd = dimsizes(dims)

   dimU = dims(nd-1)
   dimV = dims(nd-2)
   ;print("dimU = "+dimU)
   ;print("dimV = "+dimV)

; this is land (antarctica)
   varout(0,:) = varin@_FillValue

; use CICE routine to_tgrid (ice_grid.F90)
   do j = 1,dimV-1
    do i = 0,dimU-1
       ii = i-1
      if(ii .eq. -1)then
       ii = dimU-1
      end if
     val_0  = varin(j  , i)*areaU(j  , i)  ;   w----0
     val_w  = varin(j  ,ii)*areaU(j  ,ii)  ;   |    |
     val_s  = varin(j-1, i)*areaU(j-1, i)  ;   |    |
     val_sw = varin(j-1,ii)*areaU(j-1,ii)  ;   sw---s
     
     varout(j,i) = 0.25*(val_0 + val_w + val_s + val_sw)/areaT(j,i)
    end do
   end do
  assignFillValue(varin, varout)
  return(varout)
end
