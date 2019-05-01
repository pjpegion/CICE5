undef("create_outfile_name")
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
     urot@long_name=str_sub_str(urot@long_name,"X","Zonal")
     vrot@long_name=str_sub_str(urot@long_name,"Y","Meridional")
;if(1.eq.0)then
     ; delete cell_methods and standard_name attriburotes if present
;     if(isatt(urot,"cell_methods"))then
;      delete(urot@cell_methods)
;     end if
;    if(isatt(vrot,"cell_methods"))then
;      delete(vrot@standard_name)
;     end if
     if(isatt(urot,"standard_name"))then
      delete(urot@standard_name)
     end if
     if(isatt(vrot,"standard_name"))then
      delete(vrot@standard_name)
     end if
;end if

     vecs[0] = urot
     vecs[1] = vrot
     return(vecs)
end
