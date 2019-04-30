undef("create_outfile_name")
undef("regrid_and_mask")
undef(";unstagger")
undef("rotate_to_EW")

function create_outfile_name(infile:string,outres:string,outgrid:string)
    local indstr,outfile
begin
    outfile = infile
    ; find the index where the actual filename begins
     indstr = str_index_of_substr(outfile, "/", -1)
    ; insert a sub-directory name
    outfile = str_insert(outfile,outgrid+"/",indstr+1)
    ; insert a string indicating regridding resolution
    outfile = str_insert(outfile,outres,-4)
    ; find the index where the actual filename begins
     indstr = str_index_of_substr(outfile, "ocn", -1)
    ; insert an 'r' after the leading occurence of 'ocn' 
    ; to indicate regular grid
    outfile = str_insert(outfile,"r",indstr+3)
    ;print(infile)
    ;print(outfile)
    return(outfile)
end

function regrid_and_mask2d(amask:numeric,avar:numeric,wfile:string)
local amask, rgvar, avar
begin
      rgvar = new(dimsizes(amask),typeof(avar))
      amask = where(amask .gt. 0.0,amask@_FillValue,1.0)
      rgvar(0,:,:) = ESMF_regrid_with_weights(avar(0,:,:),wfile,False)

      rgvar = where(ismissing(amask),avar@_FillValue,rgvar)
    return(rgvar)
end

function regrid_and_mask3d(amask:numeric,avar:numeric,wfile:string,kdim:numeric)
local amask,avar1,rgvar,avar
begin
      rgvar = new(dimsizes(amask),typeof(avar))
      amask = where(amask .gt. 0.0,amask@_FillValue,1.0)
     do k = 0,kdim-1
       avar1 =  avar(:,k,:,:)
      rgvar(0,k,:,:) = ESMF_regrid_with_weights(avar1(0,:,:),wfile,False)
     end do
      rgvar = where(ismissing(amask),avar@_FillValue,rgvar)

    return(rgvar)
end

function unstagger(vecs:list,ufile:string,vfile:string,kdim:numeric)
local uvel,vvel,ut,vt
begin
    uvel = vecs[0]
    vvel = vecs[1]
    ; first copy Metadata
      ut = uvel
      vt = vvel
    if(dimsizes(dimsizes(ut)) .eq. 3)then
        ut(0,:,:) =   ESMF_regrid_with_weights(uvel(0,:,:),    ufile, False)
        vt(0,:,:) =   ESMF_regrid_with_weights(vvel(0,:,:),    vfile, False)
    else
     do k = 0,kdim-1
      ut(0,k,:,:) =   ESMF_regrid_with_weights(uvel(0,k,:,:),  ufile, False)
      vt(0,k,:,:) =   ESMF_regrid_with_weights(vvel(0,k,:,:),  vfile, False)
     end do
    end if
     vecs[0] = ut
     vecs[1] = vt

    return(vecs)
end

function rotate_to_EW(vecs:list,cosrot:numeric,sinrot:numeric,kdim:numeric)
local ut,vt,urot,vrot
begin
       ut = vecs[0]
       vt = vecs[1]
     ; first copy Metadata
     urot = ut
     vrot = vt
     dims = dimsizes(dimsizes(urot))
     if(dims .eq. 3)then
          urot(0,:,:) =   ut(0,:,:)*cosrot +   vt(0,:,:)*sinrot
          vrot(0,:,:) =   vt(0,:,:)*cosrot -   ut(0,:,:)*sinrot
     else
       do k = 0,kdim-1
        urot(0,k,:,:) = ut(0,k,:,:)*cosrot + vt(0,k,:,:)*sinrot
        vrot(0,k,:,:) = vt(0,k,:,:)*cosrot - ut(0,k,:,:)*sinrot
       end do
     end if

     ; change attribute to indicate these are now rotated velocities
     ut@long_name=str_sub_str(ut@long_name,"X","Zonal")
     vt@long_name=str_sub_str(vt@long_name,"Y","Meridional")
if(1.eq.0)then
     ; delete cell_methods and standard_name attributes if present
     if(isatt(ut,"cell_methods"))then
      delete(ut@cell_methods)
     end if
     if(isatt(vt,"cell_methods"))then
      delete(vt@standard_name)
     end if
     if(isatt(ut,"standard_name"))then
      delete(ut@standard_name)
     end if
     if(isatt(vt,"standard_name"))then
      delete(vt@standard_name)
     end if
end if

     vecs[0] = urot
     vecs[1] = vrot
     return(vecs)
end
