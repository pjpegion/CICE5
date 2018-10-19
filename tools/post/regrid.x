ncks -O -v lonT,latT,kmt /scratch3/NCEPDEV/marine/noscrub/Denise.Worthen/EMC_CICE/tools/grid/grid_cice_NEMS_mx025.nc tmp.nc
ncap -O -s 'land=float(kmt);land@_FillValue=0.0' tmp.nc tmp.nc
ncatted -O -a units,latT,o,c,degrees_north tmp.nc tmp.nc
ncatted -O -a units,lonT,o,c,degrees_east tmp.nc tmp.nc
mpirun -np 4 /apps/esmf/7.0.2/intel/intelmpi/bin/binO/Linux.intel.64.intelmpi.default/ESMF_RegridWeightGen -s tmp.nc --src_missingvalue land --src_coordinates lonT,latT --src_type GRIDSPEC  -d /scratch4/NCEPDEV/ocean/save/Denise.Worthen/NEMS_INPUT0.1/regrids/etopo024_oceanmask.nc --dst_missingvalue wet --dst_coordinates lon,lat --dst_type GRIDSPEC  --pole all --ignore_unmapped --ignore_degenerate -m neareststod -w cice5_tripole_tgrid_rect024.nc
rm -rf tmp.nc
