#!/bin/bash

nemsrc="/scratch4/NCEPDEV/ocean/save/Denise.Worthen/NEMS_INPUT0.1"
dstdir=$nemsrc"/regrids"
export appdir=/apps/esmf/7.0.2/intel/intelmpi/bin/binO/Linux.intel.64.intelmpi.default

export srcgrd=$WORK/EMC_CICE/tools/grid/grid_cice_NEMS_mx025.nc
export srclatlon=lonT,latT
export srcmask=kmt

# make a temporary file with missing value in mask and lat lon in degrees
# this requires that the above srcgrd file has been run in 'debug' mode so that
# the file contains the tlon and tlat
echo "ncks -O -v $srclatlon,$srcmask $srcgrd tmp.nc"
echo "ncap -O -s 'land=float(kmt);land@_FillValue=0.0' tmp.nc tmp.nc"
echo "ncatted -O -a units,latT,o,c,degrees_north tmp.nc tmp.nc"
echo "ncatted -O -a units,lonT,o,c,degrees_east tmp.nc tmp.nc"

# use the temporary file as a grid source grid
export srcgrd=tmp.nc
export srclatlon=lonT,latT
export srcmv=land

#rectilinear destination grid
export dstgrd=$nemsrc/regrids/etopo024_oceanmask.nc
export dstlatlon=lon,lat
export dstmv=wet

#export method=bilinear
export method=neareststod
export srcstring="-s $srcgrd --src_missingvalue $srcmv --src_coordinates $srclatlon --src_type GRIDSPEC "
export dststring="-d $dstgrd --dst_missingvalue $dstmv --dst_coordinates $dstlatlon --dst_type GRIDSPEC "

export weightfile="cice5_tripole_tgrid_rect024.nc"

#echo "source " $srcstring
#echo "destination " $dststring
#echo $appdir

echo "mpirun -np 4 $appdir/ESMF_RegridWeightGen $srcstring $dststring --pole all --ignore_unmapped --ignore_degenerate -m $method -w $weightfile"

echo "rm -rf tmp.nc"
