undef("create_outfile_name")
undef("regrid_and_mask")

function create_outfile_name(infile:string,outres:string,outgrid:string)
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
begin
      rgvar = new(dimsizes(amask),typeof(avar))
      amask = where(amask .gt. 0.0,amask@_FillValue,1.0)
      rgvar(0,:,:) = ESMF_regrid_with_weights(avar(0,:,:),wfile,False)

      rgvar = where(ismissing(amask),avar@_FillValue,rgvar)
    return(rgvar)
end

function regrid_and_mask3d(amask:numeric,avar:numeric,wfile:string,kdim:numeric)
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

function rotate_to_EW(vecs:list,cosrot:numeric,sinrot:numeric,kdim:numeric)
begin
    ; do ll = 0,ListCount(vecs)-1
    ;  printVarSummary(vecs[ll])
    ; end do

       ut = vecs[0]
       vt = vecs[1]
     ; first copy Metadata
     urot = ut
     vrot = vt
     dims = dimsizes(dimsizes(urot))
     if(dims .eq. 3)then
       urot(0,:,:) = ut(0,:,:)*cosrot + vt(0,:,:)*sinrot
       vrot(0,:,:) = vt(0,:,:)*cosrot - ut(0,:,:)*sinrot
     else
       do k = 0,kdim-1
        urot(0,k,:,:) = ut(0,k,:,:)*cosrot + vt(0,k,:,:)*sinrot
        vrot(0,k,:,:) = vt(0,k,:,:)*cosrot - ut(0,k,:,:)*sinrot
       end do
     end if
          ut = urot
          vt = vrot
     vecs[0] = ut
     vecs[1] = vt

     return(vecs)
end
