## Changes since 0.6 (February 13, 2014)

### Improvements
  - Implement the three-equation parameterization of the sub-shelf
    energy and salt flux balance described in *Holland D. and Jenkins
    A.*, **Modeling thermodynamic ice-ocean interactions at the base
    of an ice shelf**, Journal of Physical Oceanography, 1999
  - Allow time and space-dependent PDD temperature standard deviation.
    See [#179](https://github.com/pism/pism/issues/179).
  - Implement a parameterization of the standard deviation of the
    near-surface air temperature use by the PDD model: `Sigma = a*T + b`,
    where `T` is the mean-annual air temperature. See
    [#265](https://github.com/pism/pism/issues/265).
  - Fix RACMO mass balance in the `std-greenland` example. Impoves
    modeled ice extent. Update the corresponding User's Manual section.
  - Check the range of ice thickness during regridding and
    bootstrapping.
  - Save CTS (cold-temperate transition surface field) when `-o_size
    big` is set. See [#262](https://github.com/pism/pism/issues/262).
  - Allow statring a run with `-hydrology routing` from a result of a
    run with `-hydrology null` (set basal water amount "`bwat`" to a
    default value if not found in an input file).
  - The Temperature Index surface mass balance (PDD) model supports
    non-contiguous updates (and the `-surface ...,cache` modifier).
  - Add the `-timestep_hit_multiples X` command-line option. See
    [#254](https://github.com/pism/pism/issues/254).

### Documentation and the build system
  - Update Cray XK6 installation instructions. See
    [#240](https://github.com/pism/pism/issues/240).
  - Update Debian installation instructions.
  - Support Mac OS X 10.9.x.
  - Remove the undocumented Storglaciaren example from the User's Manual.
  - Update subglacial strength documentation.

### Bug fixes, etc
  - Add `IceModelVec2T::init_constant()` (allows using `IceModelVec2T`
    to store a field that is constant in both time and space; used by
    the PDD model).
  - Fixes in the PDD model: Fix the snow depth reset code;
    sub-intervals consistent with the interpolation scheme; fix
    initialization message (repeatable vs non-repeatable case).
  - Update lists of variables to save in PISM output files. See
    [#264](https://github.com/pism/pism/issues/264).
  - Standard out ice area and volume reporting consistent with
    saved scalar time-series.
  - Read `Href` from an input file during bootstrapping.
  - Minor improvements in model initialization code.
  - `IceGrid::kBelowHeight()` calls `MPI_Abort()`.
  - Fix a bug in `IceModelVec3Custom`.
  - Avoid re-reading periodic climate forcing data.
  - Fix a bug in the automatic vertical grid extension code.
  - Fix a bug in the code processing `-ts_times`, `-extra_times`
    command-line options: make sure that `-extra_times a:delta:b` saves
    diagnostic quantities at time b.
  - Now `-energy none` runs expect to find enthalpy, not temperature,
    in an input file.
  - Fix a bug in the `-part_redist` (redistribution of residual ice
    thickness in the sub-grid calving front position parameterization).
  - Fix a crash when reading a variable with missing "units" attribute.
  - Minor fix in `PISMRoutingHydrology`.
  - `-surface_given_period X` takes an integer X.
  - Minor fix in the discharge flux reporting code.
  - Fix a bug in `-ocean ...,cache` and `-surface ...,cache`.

## Changes since 0.5 (around May 18, 2012)

### Basal strength and basal hydrology
  - A significant change to the model physics for till: effective
	pressure in till is now computed by empirically-based (and exponential)
	formula from [`Tulaczyketal2000`] relating amount in till to effective
	pressure.
  - New mass-conserving subglacial hydrology model.

### Marine ice sheet modeling
  - Implement a parameterization of the melange back pressure as a
	part of the stress boundary condition at the ice shelf front.
  - Implement fracture density model and fracture-induced ice
    softening (see [`AlbrechtLevermann2012`]).
  - Parameterization of the sub-grid position of the grounding line and
    related improvements (see [`Feldmannetal2014`,`Gladstoneetal2012`]).
  - MISMIP3D example (see [`MISMIP3d2013`]), grounding line reversibility.
  - Fixes related to the handling of the grounding line (i.e. driving stress, etc)
  - Use an implementation of a serial two-scan connected-component
    algorithm to identify icebergs.
  - Report cumulative discharge (calving) flux as a 2D field.

### Climate inputs and ocean inputs
  - `PISMSurfaceModel` returns climatic mass balance in `[kg m-2 s-1]`
  - Remove `pclimate`; add the `-test_climate_models` option.
  - `-part_grid` extended to grounded marine termini.
  - Require time bounds in time-dependent climate forcing.
  - Improve the PDD scheme. (Keeps track of the snow depth during the year
    to choose PDD factors but does not model multi-year snow cover.)
  - Implement the 3-equation sub-shelf melt parameterization (see
    [`Hellmer1998`]).
  - Implement caching surface and ocean models `-ocean ...,cache` and
    `-surface ...,cache`. (Update boundary inputs every X years but
    include interannual variability.)
  - Remove `EISMINT-Greenland`.

### Inverse modeling
  - Inverse modeling tools *are* a part of the release now.

### Energy
  - New bootstrapping heuristic filling ice temperature at depth from surface
    mass balance (if available).
  - Allow "regridding" enthalpy from files containing ice temperature
    and liquid water fraction.
  - Allow cold flow laws in the enthalpy mode and GPBLD in the cold
    mode (same as Paterson-Budd).
  - Corrected basal boundary condition in the enthalpy system.

### Usability
  - Implement poor man's parallel I/O, with compression
  - Ensure "continuity" in time of reported cumulative diagnostic fields.
  - Let the user precisely specify the dates corresponding to the run (`-time_file`).
  - Many more diagnostic quantities; use `-list_diagnostics` to see.
  - PISM keeps track of options and parameters that were set but were not used.
  - Use projection info to compute latitude and longitude bounds.
	(Reduces post-processing needed to work with PISM's output.)

### Miscellaneous
  - Improve mass conservation and mass transport.
  - New validation example based on laboratory experiment with xanthan gum.
    Uses millimeter grid spacing and shows configurability of PISM.
  - Implement a "constant 2D velocity" stress balance object (for testing and
    prescribing sliding)

### Under the hood
  - Improve the build system.
  - Switch to PETSC 3.3 or 3.4, stop supporting 3.2.
  - Switch to [UDUNITS-2](http://www.unidata.ucar.edu/software/udunits/).
  - Require [FFTW-3](http://www.fftw.org) to build PISM.
  - Use [CalCalcs](http://meteora.ucsd.edu/~pierce/calcalcs/) for proper calendar support.
  - Remove the Debian [meta-] package.
  - Clean up command-line options selecting sub-models (calving,
    stress balance, energy, basal yield stress).
  - Better `SSAFD` recovery logic; save model state on failure.
  - Use non-zero initial guess in the SSAFD KSP solver.
  - Improved basal yield stress code.

## Changes from 0.4 (around June 1, 2011) to 0.5

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
