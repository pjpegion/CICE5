#! /bin/csh -f

set CDATE="20160101"
set dirname="/scratch4/NCEPDEV/nems/noscrub/Bin.Li/benchmark_test/$CDATE/COMFV3/c384_test/gfs.$CDATE/00/ICE/"

set hemi="NH"
set outname="ice."$hemi
ncl dirname=\"{$dirname}\"  'varnames = (/"aice", "hi", "hs", "albsni"/)' hemi=\"{$hemi}\" outname=\"{$outname}\" benchplots_scalar.ncl
set outname="melt.pond."$hemi
ncl dirname=\"{$dirname}\"  'varnames = (/"meltt", "meltb", "apond", "hpond"/)' hemi=\"{$hemi}\" outname=\"{$outname}\" benchplots_scalar.ncl

set hemi="SH"
set outname="ice."$hemi
ncl dirname=\"{$dirname}\"  'varnames = (/"aice", "hi", "hs", "albsni"/)' hemi=\"{$hemi}\" outname=\"{$outname}\" benchplots_scalar.ncl
set outname="melt.pond."$hemi
ncl dirname=\"{$dirname}\"  'varnames = (/"meltt", "meltb", "apond", "hpond"/)' hemi=\"{$hemi}\" outname=\"{$outname}\" benchplots_scalar.ncl
