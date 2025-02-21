#!/usr/bin/env python3
#
# Copyright (C) 2025 Andy Aschwanden
#
# This file is part of PISM.
#
# PISM is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# PISM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License
# along with PISM; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

# PISM-Storglaciaren 3D example
#
# Here we apply PISM to a small mountain glacier in northern Sweden,
# Storglaciaren, using a grid refinement path for the shallow-ice
# approximation (SIA), the hybrid approximation (HY), and the Higher-Order
# approximation to the Stokes Equations. We make use of PISM's grid-description
# capabilities, that is, we define the modeling domain via a CDL file which is
# converted to a netCDF file and prescribe the horizontal grid resolution in meters
# instead of number of grid points.
# We choose ice flow parameters such that the higher-order model predicts basal and
# surface speeds in broad agreement with observations while the shallow ice and
# hybrid model predict speeds higher than observations.


# prefix to pism (not to executables)
if [ -n "${PISM_PREFIX:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME     PISM_PREFIX = $PISM_PREFIX  (already set)"
else
  PISM_PREFIX="$HOME/pism"    # just a guess
  echo "$SCRIPTNAME     PISM_PREFIX = $PISM_PREFIX"
fi

# Run the preprocssing script that generates the initial state file
./preprocess.sh

PISM=$PISM_PREFIX/bin/pism

bootfile=pism_storglaciaren_3d.nc
gridfile=grid.nc

# setting up the SIA, Hybrid, and Higher-Order Solver
sia="-skip -skip_max 250 -stress_balance sia -stress_balance.sia.enhancement_factor 0.3"
hy="-skip -skip_max 250 -stress_balance ssa+sia -stress_balance.sia.enhancement_factor 0.3"
ho="-stress_balance blatter -stress_balance.blatter.enhancement_factor 0.3 -stress_balance.blatter.coarsening_factor 4 -blatter_Mz 17 -bp_ksp_type gmres -bp_pc_type mg -bp_pc_mg_levels 3 -bp_mg_levels_ksp_type richardson -bp_mg_levels_pc_type sor -bp_mg_coarse_ksp_type preonly -bp_mg_coarse_pc_type lu -bp_snes_monitor_ratio -bp_ksp_monitor -bp_ksp_view_singularvalues -bp_snes_ksp_ew 1 -bp_snes_ksp_ew_version 3 -time_stepping.adaptive_ratio 250"

# basal conditions for Hybrid and Higher-Order
basal="-basal_resistance.pseudo_plastic.enabled yes  -basal_resistance.pseudo_plastic.q 0.75 -basal_resistance.pseudo_plastic.u_threshold 1.0 -basal_yield_stress.mohr_coulomb.till_phi_default 30 -stress_balance.sia.bed_smoother.range 10"

# Climate forcing
climate="-surface given,forcing -surface_given_file $bootfile -surface.force_to_thickness.file $bootfile"

regridvars="litho_temp,enthalpy,age,tillwat,bmelt,ice_area_specific_volume,thk"

    
declare -a sbs=("$sia" "$hy" "$ho")
declare -a prefixs=("sia" "hy" "ho")

n=${#sbs[@]}
for (( i=1; i<${n}+1; i++ )); do
    sb=${sbs[$i-1]}
    prefix=${prefixs[$i-1]}
    
    N=2
    res=80m
    y=1000
    outfile=sg_${prefix}_${res}_${y}a.nc
    exfile=ex_$outfile
    tsfile=ts_$outfile
    grid="-grid.dx $res -grid.dy $res -Mbz 1 -Lz 500 -Mz 51 -z_spacing equal -grid.file $gridfile"
    cmd="mpirun -np $N $PISM -config_override psg_config.nc -bootstrap -i $bootfile $climate -y $y  $grid  -extra_times 10 -extra_vars usurf,taud_mag,topg,tauc,taud,beta,mask,thk,temppabase,tempicethk_basal,velbase_mag,velsurf_mag,climatic_mass_balance,tillphi -output.extra.stop_missing no -extra_file $exfile -ts_times yearly -ts_file $tsfile -o $outfile $basal $sb"
    $cmd

    N=4
    res=40m
    y=100
    infile=$outfile
    outfile=sg_${prefix}_${res}_${y}a.nc
    exfile=ex_$outfile
    tsfile=ts_$outfile
    grid="-grid.dx $res -grid.dy $res -Mbz 1 -Lz 500 -Mz 51 -z_spacing equal -grid.file $gridfile"
    cmd="mpirun -np $N $PISM -config_override psg_config.nc -bootstrap -i $bootfile -input.regrid.file $infile -input.regrid.vars $regridvars  $climate -y $y $grid  -extra_times 10 -extra_vars usurf,taud_mag,topg,tauc,taud,beta,mask,thk,temppabase,tempicethk_basal,velbase_mag,velsurf_mag,climatic_mass_balance,tillphi -output.extra.stop_missing no -extra_file $exfile -ts_times yearly -ts_file $tsfile -o $outfile $basal $sb"
    $cmd

    N=8
    res=20m
    y=10
    infile=$outfile
    outfile=sg_${prefix}_${res}_${y}a.nc
    exfile=ex_$outfile
    tsfile=ts_$outfile
    grid="-grid.dx $res -grid.dy $res -Mbz 1 -Lz 500 -Mz 51 -z_spacing equal -grid.file $gridfile"
    cmd="mpirun -np $N $PISM -config_override psg_config.nc -bootstrap -i $bootfile -input.regrid.file $infile -input.regrid.vars $regridvars  $climate -y $y $grid  -extra_times 1 -extra_vars usurf,taud_mag,topg,tauc,taud,beta,mask,thk,temppabase,tempicethk_basal,velbase_mag,velsurf_mag,climatic_mass_balance,tillphi -output.extra.stop_missing no -extra_file $exfile -ts_times yearly -ts_file $tsfile -o $outfile $basal $sb"
    $cmd
done

python plot_results.py sg_*_80m_1000a.nc
python plot_results.py sg_*_40m_100a.nc
python plot_results.py sg_*_20m_10a.nc
