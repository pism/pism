## Changes since 0.4 (around June 1, 2011)

  - Switch to PETSc 3.2
  - Move to github.com
  - Add regional modeling tools: pismo, Python tools
  - Improve the IceFlowLaw class by moving physical constants out of it
  - Implement `-o_format [netcdf4_parallel, pnetcdf]`
  - Implement `-o_order [xyz, yxz, zyx]`
  - Implement `-e_age_coupling`
  - Improved the implementation of PISM's basal yield stress model
  - Separated flow laws used in SIA and SSA code. New command-line
    options: `-sia_e`, `-ssa_e`. New config. parameters:
    `sia_enhancement_factor`, `ssa_enhancement_factor`
  - Fixed `-topg_to_phi`: now it takes 4 numbers (no special value in the ocean)
  - `-pdd_annualized` implemented, and recommended for degree-day scheme usage
  - Flush time-series when `-extra_files` are written
  - Improved documentation (mostly climate forcing code)
  - Improved climate forcing code and its interface (in particular:
    consistent command-line options)
  - Fixed numerous bugs
  - Improved regridding code
  - Implemented temperature-dependent conductivity and specific heat capacity of ice
  - Model time is in seconds
  - PISM uses time bounds, both in reporting code and code using forcing data
  - Reported rates of change are computed as average rates over
    reporting intervals, by differencing cumulative quantities
  - Reporting: daily, monthly, yearly
  - Added support for the Gregorian calendar
  - PISM uses Proj.4 to compute latitudes, longitudes, and cell areas
  - Implemented `-tauc_to_phi`
  - Updated and improved modeling examples
  - Removed the pgrn executable; all whole-Greenland examples use pismr
  - Removed the EISMINT-Ross example and the pross executable
  - Inverse modeling tools (not part of the release)

## Changes since 0.5 (around May 18, 2012)

### Basal strength and basal hydrology
  - A significant change to the model physics for till: effective
	pressure in till is now computed by empirical (and exponential)
	formula from [`Tulaczyketal2000] relating amount in till to effective
	pressure.
  - New subglacial hydrology model

### Marine ice sheet modeling
  - Implement a parameterization of the melange back pressure as a
	part of the stress boundary condition at the ice shelf front
  - Implement fracture density model and fracture-induced ice
    softening (see [`AlbrechtLevermann2012`])
  - Work on the sub-grid grounding line position approximation and
    related improvements (see [`Feldmannetal2014`])
  - MISMIP3D (see [`MISMIP3d2013`]), grounding line reversibility
  - Fixes related to the handling of the grounding line (i.e. driving stress, etc)
  - Use an implementation of a serial two-scan connected-component
    algorithm to identify icebergs.
  - Report cumulative discharge (calving) flux as a 2D field

### Climate inputs and ocean inputs
  - `PISMSurfaceModel` returns climatic mass balance in `[kg m-2 s-1]`
  - Remove `pclimate`; add the `-test_climate_models` option
  - `-part_grid` extended to grounded marine termini
  - Fix and clean up the parameterization of the sub-grid position of
    the grounding line (see [`Gladstoneetal2012`])
  - Require time bounds in time-depencent climate forcing.
  - Improve the PDD scheme. (Keeps track of the snow depth to during
    the year to choose PDD factors but does not model multi-year snow
    cover.)
  - Implement the 3-equation sub-shelf melt parameterization (see
    Hellmer1998)
  - Implement caching surface and ocean models `-ocean ...,cache` and
    `-surface ...,cache`. (Update boundary inputs every X years but
    include interannual variability.)
  - Remove `EISMINT-Greenland`.

### Inverse modeling
  - Inverse modeling tools *are* a part of the release now.

### Energy
  - Implement the new bootstrapping heuristic (filling ice temperature at depth)
  - Allow "regridding" enthalpy from files containing ice temperature
    and liquid water fraction
  - Allow cold flow laws in the enthalpy mode and GPBLD in the cold
    mode (same as PB).
  - Better understood basal B.C. in the enthalpy system. Clean up the enthalpy code.

### Usability
  - Implement poor man's parallel I/O, with compression
  - Ensure "continuity" in time of reported cumulative diagnostic fields
  - Let the user precisely specify the dates corresponding to the run (`-time_file`)
  - many more diagnostic quantities; `-list_diagnostics`
  - PISM keeps track of options and parameters that were set but were not used
  - Use projection info to compute latitude and longitude bounds
	(reduces post-processing needed to work with PISM's output)

### Miscellaneous
  - Improve mass conservation and mass transport
  - New validation example; shows that PISM can be used at different
	 spatial scales (gum experiment vs. Antarctica)
  - Implement a "constant 2D velocity" stress balance object (can be
    used for testing and prescribing sliding)

### Under the hood
  - Improve the build system
  - Switch to PETSC 3.3 or 3.4, stop supporting 3.2
  - Switch to [UDUNITS-2](http://www.unidata.ucar.edu/software/udunits/)
  - Require [FFTW-3](http://www.fftw.org) to build PISM
  - Use [CalCalcs](http://meteora.ucsd.edu/~pierce/calcalcs/) for proper calendar support
  - Remove the Debian [meta-] package
  - Clean up command-line options selecting sub-models (calving,
    stress balance, energy, basal yield stress)
  - Better `SSAFD` recovery logic; save model state on failure
  - Use non-zero initial guess in the SSAFD KSP solver.
  - Improved basal yield stress code
