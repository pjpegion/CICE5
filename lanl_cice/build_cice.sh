#! /bin/sh -l

echo ""
echo "=================================================================="
echo "=   build LANL CICE for NEMS                                     ="
echo "=   SYNTAX:                                                      ="
echo "=     build_cice.sh                                              ="
echo "=                                                                ="
echo "=================================================================="


cwd=$(pwd)

source ~/bin/loadmodule.sh
echo "Netcdf library: $NETCDF"

echo "Changing direction to $HOME/github/mom/exp"
echo ""
cd $HOME/github/lanl_cice
./comp_ice
ar -r /home/Fei.Liu/noscrub/lanl_cice/libcice.a CICE_FinalMod.o CICE_InitMod.o CICE_RunMod.o ice_aerosol.o ice_age.o ice_algae.o ice_atmo.o ice_blocks.o ice_boundary.o ice_brine.o ice_broadcast.o ice_calendar.o ice_communicate.o ice_constants.o ice_diagnostics.o ice_distribution.o ice_domain.o ice_domain_size.o ice_dyn_eap.o ice_dyn_evp.o ice_dyn_shared.o ice_exit.o ice_fileunits.o ice_firstyear.o ice_flux.o ice_forcing.o ice_gather_scatter.o ice_global_reductions.o ice_grid.o ice_history.o ice_history_bgc.o ice_history_drag.o ice_history_mechred.o ice_history_pond.o ice_history_shared.o ice_history_write.o ice_init.o ice_itd.o ice_kinds_mod.o ice_lvl.o ice_mechred.o ice_meltpond_cesm.o ice_meltpond_lvl.o ice_meltpond_topo.o ice_ocean.o ice_orbital.o ice_read_write.o ice_restart.o ice_restart_driver.o ice_restart_shared.o ice_restoring.o ice_shortwave.o ice_spacecurve.o ice_state.o ice_step_mod.o ice_therm_0layer.o ice_therm_bl99.o ice_therm_itd.o ice_therm_mushy.o ice_therm_shared.o ice_therm_vertical.o ice_timers.o ice_transport_driver.o ice_transport_remap.o ice_zbgc.o ice_zbgc_shared.o shr_orb_mod.o
cd $cwd

module list
