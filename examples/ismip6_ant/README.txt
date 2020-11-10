!---------------------------------------------------------------------------------
! Test for ISMIP6 oceanic parameterisation from Jourdain et al. (2020).
!
! Those scripts set up a short bootstrapping run and a short no-mass run to
! test the ISMIP6 oceanic models for Antarctica.
!
!---------------------------------------------------------------------------------

User needs to adapt the script to their own environment:
- set_environment.sh: local user environment
- run_paleo.sh: location of input files and PISM configuration file need to be adjusted as well
- run_pismr.sh: location of directories to be modified

ISMIP6 oceanic models are:
----------------------------
- ocean.ismip6: local Parameterisation
- ocean.ismip6nl: non-local Parameterisation

They can be set in: set_physics.sh

Input files:
---------------
- atmospheric and topographic conditions: bedmap2_martos_MAR_monthly_pism30km.nc
- oceanic temperature, salinity and basins: schmidtko_pism30km_means_ismip6.nc

The oceanic file can be used by both ISMIP6 models,
the one not using drainage basin (local parameterisation) ignores the basin_mask variable.


RUNNING:
---------------

 sh run_pismr.sh

 to run all the scripts.
