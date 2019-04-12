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

function regrid_and_mask(amask:numeric,avar:numeric,wfile:string)
begin
      rgmask = ESMF_regrid_with_weights(amask,wfile,False)
      rtmask = where(rgmask .gt. 0.0,avar@_FillValue,1.0)

      rgvar = ESMF_regrid_with_weights(avar,wfile,False)
      rgvar = where(ismissing(rtmask),avar@_FillValue,rgvar)
    return(rgvar)
end


