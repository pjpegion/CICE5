; This is a copy of the curvilinear_to_SCRIP script provided 
; inside of EMSF_regridding.ncl as of NCL v 6.5.0, 30-April-2019
; which does not require the grid input values to be re-ordered
; 
;======================================================================
; This function receives a 2D curvilinear grid and mask and stores
; it in FName NetCDF file based on SCRIP standard.
; lat2d and lon2d must have (nlat,nlot) dimension.
; Opt controls the behavior of the function
; current attribute in Opt are:
; (1) Overwrite [Logical] if True, and if the out file already exists
;     it will erase the file.
; (2) ForceOverwrite [Logical] If set to True, the user is not asked
;     for removing an existing file. If set to false the user permission
;     is required to remove an existing file. This is ineffective if
;     Overwrite is set to False.
; (3) If GridCornerLat and GridCornerLon are set, then these will be used
;     instead of calcuation the corner points for the lat/lon cells.
;     The corner points are needed for the "conserve" method.
;======================================================================
undef("my_curvilinear_to_SCRIP")
procedure my_curvilinear_to_SCRIP(FName[1]:string,lat2d[*][*]:numeric,\
                               lon2d[*][*]:numeric,Opt[1]:logical)
local latlon_dims, nlat, nlon, fid, FileAtt, grid_siz, grid_corners, \
grid_rank, FDimNames, FDimSizes, FDimUnlim, DummyAtt1, DummyAtt2, \
GridCornerLat, GridCornerLon, grid_corner_lat, grid_corner_lon, grid_mask, \
lat_type, lon_type, DEBUG, nc_file_type
begin
;---Check for options
    PrintTimings = isatt_logical_true(Opt,"PrintTimings")
    if(PrintTimings) then
      start_time = get_start_time()
    end if

    DEBUG = isatt_logical_true(Opt,"Debug")
    TESTIT = isatt_logical_true(Opt,"Testit")

;---Check if the file already exists
    check_for_file(FName,Opt)

;---Do we need to create a large file?
    if (.not.isatt_logical_false(Opt,"LargeFile")) then
      setfileoption("nc","Format","LargeFile") 
    end if

;---Do we need to create a netcdf-4 file?
    if (isatt(Opt,"NetCDFType").and.\
        str_lower(Opt@NetCDFType).eq."netcdf4") then
      setfileoption("nc","Format","netcdf4") 
      nc_file_type = "netcdf4"
    else
      nc_file_type = "netcdf3"
    end if

    if ( any(dimsizes(lat2d).ne.dimsizes(lon2d)) ) then
        print("curvilinear_to_SCRIP: latitude and longitude must have the same number of elements.")
        exit
    end if

    latlon_dims = dimsizes(lat2d)
    nlat        = latlon_dims(0)
    nlon        = latlon_dims(1)

    if(Opt.and.isatt(Opt,"Title")) then
      FTitle = Opt@Title
    else
      FTitle = "curvilinear_to_SCRIP (" + nlat + "," + nlon + ")"
    end if

;---Was a mask provided?
    grid_mask_name = get_mask_name(Opt)
    if(grid_mask_name.ne."") then
      if(.not.all(dimsizes(Opt@$grid_mask_name$).eq.latlon_dims)) then
        print("curvilinear_to_SCRIP: Opt@" + grid_mask_name + \
              " is not the correct dimensionality")
        exit
      else
        grid_mask = Opt@$grid_mask_name$
      end if
    else
;---No masking
      grid_mask = new(latlon_dims, "integer","No_FillValue")
      grid_mask = 1
    end if

;---Create the file
    fid = addfile(FName,"c")
    setfileoption(fid,"DefineMode",True)

;---Define the file attributes
    FileAtt              = True
    FileAtt@title        = FTitle
    FileAtt@Conventions  = "SCRIP"
    FileAtt@Createdby    = "ESMF_regridding.ncl"
    FileAtt@date_created = systemfunc("date")
    fileattdef(fid,FileAtt)
    
;---Define the SCRIP dimensions
    grid_size    = nlat*nlon   ; This is number of data points (grid nodes)
    grid_corners = 4	
    grid_rank    = 2
	
    FDimNames = (/ "grid_size","grid_corners","grid_rank" /)
    FDimSizes = (/ grid_size,grid_corners,grid_rank /)
    FDimUnlim = (/ False,False,False /)
	
    filedimdef(fid,FDimNames,FDimSizes,FDimUnlim)
	
;---Define Variables
    filevardef(fid,"grid_dims","integer","grid_rank")
    filevardef(fid,"grid_center_lat","double","grid_size")
    filevardef(fid,"grid_center_lon","double","grid_size")
    filevardef(fid,"grid_imask","integer","grid_size")
    filevardef(fid,"grid_corner_lat","double",(/ "grid_size", "grid_corners" /) )
    filevardef(fid,"grid_corner_lon","double",(/ "grid_size", "grid_corners" /) )

;---Define the variables unit attribute
    DummyAtt1       = 0
    DummyAtt2       = 0
;
; The conversion to character arrays for NetCDF-4 files is a 
; work-around.  ESMF_RegridWeightGen uses the F90 interface to 
; read NetCDF-4 files, which doesn't yet have support for 
; NetCDF-4 "string" types.
; 
    if(nc_file_type.eq."netcdf3") then
      DummyAtt1@units = "degrees"
      DummyAtt2@units = "unitless"
    else
      DummyAtt1@units = tochar("degrees")
      DummyAtt2@units = tochar("unitless")
    end if
    filevarattdef(fid,"grid_center_lat",DummyAtt1)
    filevarattdef(fid,"grid_center_lon",DummyAtt1)
    filevarattdef(fid,"grid_imask",DummyAtt2)
    filevarattdef(fid,"grid_corner_lat",DummyAtt1)
    filevarattdef(fid,"grid_corner_lon",DummyAtt1)
    delete(DummyAtt1)
    delete(DummyAtt2)
    
;---Prepare the file to store the values
    setfileoption(fid,"DefineMode",False)
    
;---Storing Grid Dims
    fid->grid_dims = (/ nlon, nlat /)   ; SCRIP is FORTRAN-based.
                                        ; (nlat,nlon) in NCL is equivalent to
                                        ; (nlon,nlat) in FORTRAN
    
;---Store Cell Center Lat/Lon
    fid->grid_center_lat = (/ndtooned(lat2d)/)
    fid->grid_center_lon = (/ndtooned(lon2d)/)
	
;---Store Cell Masks
    if (grid_size.ne.dimsizes(ndtooned(grid_mask))) then
        print("curvilinear_to_SCRIP: Mask array is not the appropriate size.")
        exit
    else
        fid->grid_imask=(/ tointeger(ndtooned(grid_mask)) /)
    end if
	
;---Get the grid lat/lon corners
    if (isatt(Opt,"GridCornerLat").and.isatt(Opt,"GridCornerLon")) then
        if(DEBUG) then
          print("curvilinear_to_SCRIP: using grid corners provided by user...")
        end if
        if(.not.any(typeof(Opt@GridCornerLat).eq.(/"float","double"/))) then
          GridCornerLat = tofloat(Opt@GridCornerLat)
        else
          GridCornerLat = Opt@GridCornerLat
         ; dims = dimsizes(GridCornerLat)
         ; print("GridCornerLat dims  "+dims)
         ; delete(dims)
         ; print(GridCornerLat(:,0,0))
        end if
        if(.not.any(typeof(Opt@GridCornerLon).eq.(/"float","double"/))) then
          GridCornerLon = tofloat(Opt@GridCornerLon)
        else
          GridCornerLon = Opt@GridCornerLon
        end if
        ;grid_corner_lat = reshape( GridCornerLat,(/ grid_size, grid_corners /))
        ;grid_corner_lon = reshape( GridCornerLon,(/ grid_size, grid_corners /))
      
        ;dims = dimsizes(grid_corner_lat)
        ;print("grid_corner_lat dims  "+dims(0)+" "+dims(1))
; if reshape is doing what it is supposed to do, this should print same 4 values
        ;print(grid_corner_lat(0,:))
        ;delete(dims)

; explicitly do what reshape is trying to do  
      grid_corner_lat = new( (/ grid_size, grid_corners/),typeof(GridCornerLat))
      grid_corner_lon = new( (/ grid_size, grid_corners/),typeof(GridCornerLon))
      do n = 0,grid_corners-1
       grid_corner_lat(:,n) = ndtooned(GridCornerLat(n,:,:))
       grid_corner_lon(:,n) = ndtooned(GridCornerLon(n,:,:))
      end do
    else
;---Estimate the grid cell corners; it's better if the user provides them.
        if(DEBUG) then
          print("curvilinear_to_SCRIP: calculating grid corners...")
        end if
        lat_type = typeof(lat2d)
        lon_type = typeof(lon2d)
        if(.not.any(lat_type.eq.(/"float","double"/))) then
          lat_type = "float"
        end if
        if(.not.any(lon_type.eq.(/"float","double"/))) then
          lon_type = "float"
        end if
        grid_corner_lat = new( (/ grid_size, grid_corners /), lat_type)
        grid_corner_lon = new( (/ grid_size, grid_corners /), lon_type)
        if(any(fabs(lat2d).eq.90)) then
          if(DEBUG) then
            print("curvilinear_to_SCRIP: one or more lat values are at the")
            print("       poles, so calculating grid corners using")
            print("       calc_SCRIP_corners_boundaries...")
          end if
          if(TESTIT) then
            calc_SCRIP_corners_boundaries_2(lat2d,grid_corner_lat,"lat",Opt)
            calc_SCRIP_corners_boundaries_2(lon2d,grid_corner_lon,"lon",Opt)
          else 
            calc_SCRIP_corners_boundaries(lat2d,lon2d,grid_corner_lat,grid_corner_lon,Opt)
          end if
        else
          if(DEBUG) then
            print("curvilinear_to_SCRIP: no lat values are at the poles, so")
            print("       calculating grid corners using")
            print("       calc_SCRIP_corners_noboundaries...")
          end if
          calc_SCRIP_corners_noboundaries(lat2d,lon2d,grid_corner_lat,grid_corner_lon,Opt)
         if(any(fabs(grid_corner_lat).gt.90)) then
           if(DEBUG) then
             print("curvilinear_to_SCRIP: calc_SCRIP_corners_noboundaries")
             print("       produced out-of-range latitude values.")
             print("       Trying calc_SCRIP_corners_boundaries...")
           end if
            if(TESTIT) then
              calc_SCRIP_corners_boundaries_2(lat2d,grid_corner_lat,"lat",Opt)
              calc_SCRIP_corners_boundaries_2(lon2d,grid_corner_lon,"lon",Opt)
            else 
              calc_SCRIP_corners_boundaries(lat2d,lon2d,grid_corner_lat,grid_corner_lon,Opt)
            end if
          end if
        end if
    end if
   
;---Store the cell corners
    if (isatt(grid_corner_lat,"_FillValue")) then
        delete(grid_corner_lat@_FillValue)
    end if
    if (isatt(grid_corner_lon,"_FillValue")) then
        delete(grid_corner_lon@_FillValue)
    end if
    fid->grid_corner_lat = (/ todouble(grid_corner_lat) /)
    fid->grid_corner_lon = (/ todouble(grid_corner_lon) /)

    if(PrintTimings) then
      print_elapsed_esmf_time(start_time,"curvilinear_to_SCRIP")
    end if
end     ; of curvilinear_to_SCRIP(...)
