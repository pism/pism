.. default-role:: literal


Changes since v2.1
==================

- Add some automatic testing on macOS.
- Use pkg-config to look for all the dependencies that support it.
- Add `VERSION` to ensure that tarballs with PISM's sources (e.g. archived by Zenodo)
  contain appropriate version info.
- Update some examples in `examples/marine`.
- Fix an bug reported by Ken Mankoff: scripts in `examples/antarctica/` required PISM
  built with PROJ.

Changes since v2.0
==================

- Support 2D precipitation offsets in `-atmosphere ...,delta_P`. If the input file set
  using `atmosphere.delta_P.file` contains a scalar time series `delta_P`, use that as a
  time-dependent constant-in-space forcing. If the input file contains a 2D variable
  `delta_P`, use that as a time-and-space-dependent forcing.
- Support 2D air temperature offsets in `-atmosphere ...,delta_T`. If the input file set
  using `atmosphere.delta_T.file` contains a scalar time series `delta_T`, use that as a
  time-dependent constant-in-space forcing. If the input file contains a 2D variable
  `delta_T`, use that as a time-and-space-dependent forcing.
- Refactor utility classes used to store 2D and 3D arrays.
- Remove a misguided energy conservation attempt that turned out to be harmful
  (occasionally).
- Fix a bug in the code reading periodic time-dependent forcing.
- Update pre-processing scripts in `examples/std-greenland`.
- Fix the scaling of the `uplift` diagnostic in `-bed_def given`.
- Fix GSL-related build issues (unable to find GSL when it is installed in a non-standard
  location).
- Fix documentation of `...till_effective_fraction_overburden`.
- Use `realpath()` to resolve relative file names. Now configuration parameters ending in
  `.file`, when saved to output files and in PISM output to `stdout`, contain *absolute*
  file names. This will make it easier to reproduce runs based on an output file.
- Support checkpointing the HTCondor way (see commit 3740c41df).
- Stop with an error message if a NetCDF variable in an input file contains not-a-number
  or infinity.
- Use `-list_diagnostics all` to print the list of all diagnostics, `-list_diagnostics
  spatial` for 2D and 3D variables, and `-list_diagnostics scalar` for scalar time series.
- Support piecewise-constant temporal interpolation of near-surface air temperatures in
  `-atmosphere given`: set `atmosphere.given.air_temperature_interpolation` to
  `piecewise_constant`.
- Extrapolate sliding velocities computed by the SSAFD solver to improve the initial guess
  used when the ice front advances (set
  :config:`stress_balance.ssa.fd.extrapolate_at_margins` to `false` to disable).

Changes since v1.2
==================

- The three-equation ocean model `-ocean th` uses constant salinity (see
  `constants.sea_water.salinity`) if `salinity_ocean` is not present in the forcing file.
- `fill_missing_petsc.py` uses homogeneous Neumann BC at domain boundaries.
- Support 2D precipitation scaling in `-atmosphere ...,frac_P`. If the input file set
  using `atmosphere.frac_P.file` contains a scalar time series `frac_P`, use that as a
  time-dependent constant-in-space forcing. If the input file contains a 2D variable
  `frac_P`, use that as a time-and-space-dependent forcing.
- Add a new `output.format` value: `netcdf4_serial` and `output.compression_level`. Use
  `-o_format netcdf4_serial -output.compression_level N` (`N` between 1 and 9) to write
  compressed NetCDF from rank 0.
- Support writing compressed NetCDF in parallel with NetCDF 4.7.4 or newer and HDF5 1.10.3
  or newer. Set `output.compression_level` to enable compression.
- Stop with an error message if some values of a variable read from a file match the
  `_FillValue` attribute (PISM expects input files to contain data at all grid points
  within the domain).
- Add optional arguments `time_units` and `calendar` to `PISM.util.prepare_output()` in
  the Python bindings.
- Add surface elevation smoothing to the orographic precipitation model. High-frequency
  modes in the surface elevation that can develop in runs with evolving ice geometry
  (consider grounded ice margins) may cause oscillations in the computed precipitation
  field (probably due to the Gibbs phenomenon). These oscillations may result in an even
  rougher topography, triggering a feedback loop polluting model results. Set
  `atmosphere.orographic_precipitation.smoothing_standard_deviation` (in meters)
  to smooth the ice surface elevation to reduce this effect.
- Add `sea_level.constant.value`. This sets the default sea level elevation used with
  `-sea_level constant`.
- Remove `ocean.always_grounded`. Set `sea_level.constant.value` to a large negative value
  to ensure that all ice is grounded.
- Remove `ocean.melange_back_pressure_fraction`: it is no longer needed.
- Add a new ocean modifier: `-ocean ...,delta_MBP`. This component reads scalar
  time-dependent melange pressure offsets (units: Pa) and uses them in the calving front
  boundary condition for the SSA.
- The new `-bed_def given` class reads in a prescribed bed deformation history from a file
  (e.g. from a solid-Earth model) relative to a (high-resolution) reference topography,
  indicated by configuration parameter `bed_deformation.given.file` and
  `bed_deformation.given.reference_file`, respectively.
- Implemented regularized Coulomb sliding as in Zoet & Iverson, 2020, A slip law for
  glaciers on deformable beds, equation 3.
- Implement a FEM solver for the first order approximation of the Stokes equations due to
  Blatter (1995). This solver supports multigrid preconditioners (see Brown et al 2013)
  and includes 5 verification test based on manufactured solutions.
- Implement experiments A,B,C,D,E from the ISMIP-HOM intercomparison.
- Adjust PICO ocean input average across covered basins, in which the ice shelf has
  in fact a connection to the ocean. Large ice shelves, that cover across two basins,
  that do not share an ocean boundary, are split into two separate ice shelf instances
- Implement scaling of calving rates using a time-dependent factor. Set
  `calving.rate_scaling.file` to the name of the file containing `frac_calving_rate`
  (units: "1").
- Add the new command-line option `-refinement_factor N`. Use this to select a regional
  modeling domain using `-x_range ... -y_range ...`, then use a grid that is `N` times
  finer.
- Fix a bug in the code managing time step restrictions (this affected the last time step
  of runs using `-skip` and runs with `-skip` in which `-max_dt` is active).
- Add support for automatic unit conversion in command-line options. If an option argument
  is a number PISM assumes that it uses PISM's internal units. If it is a number followed
  by a units string recognized by UDUNITS it is automatically converter to PISM's internal
  units. For example, the following are equivalent: `-Lz 1000`, `-Lz 1000m`, `-Lz 1km`.
- Command-line options `-y`, `-ys`, `-ye`, `-max_dt` and corresponding configuration
  parameters use units of `365 days` instead of `years`. The latter has the meaning of the
  mean tropical year, i.e. the constant used to convert from `1/s` to `1/year`. Use `-y
  1000years`, etc to reproduce the old behavior.
- Add a new parameter: `time_stepping.resolution`. PISM rounds time step lengths *down* to
  a multiple of this number (default: 1 second). This reduces the influence of rounding
  errors on time step lengths. See `issue 407`_.
- Remove the `pisms` executable. Run `pismr -eisII X` to run EISMINT-II test `X`.
- Ice thickness threshold read in from `calving.thickness_calving.file` can be both space-
  and time-dependent.
- Now PISM stops with an error message if time-dependent forcing data read from a file do
  not span the whole length of a simulation. Set `input.forcing.time_extrapolation` to
  "true" to disable this check.
- Remove the configuration parameter `input.forcing.evaluations_per_year`. Now
  the code evaluates *exact* values of time averages of time-dependent forcing inputs.
- Major improvement in the handling of time-dependent forcing. A file containing periodic
  forcing has to contain *exactly* one period. The start and the length of the period are
  derived from time bounds read from this file. This makes it easier to use periodic
  forcing and adds supports for arbitrary period lengths. See the manual section about
  time-dependent inputs.
- All time-dependent forcing files have to contain time bounds.
- Now PISM always respects the reference date in input files.
- Rename NetCDF variables `bc_mask` to `vel_bc_mask` and `u_ssa_bc` and `v_ssa_bc` to
  `u_bc` and `v_bc`.
- Add a new NetCDF variable `thk_bc_mask` prescribing locations where the ice thickness is
  kept fixed. This mask is combined with `vel_bc_mask`: we keep ice thickness fixed at all
  the locations where the sliding (usually SSA) velocity is fixed.
- Implement the fracture density growth parameterization due to Borstad et al
  (equation 4 in http://doi.org/10.1002/2015GL067365). Code contributed by T. Albrecht).
- Assume that in the "ocean" areas the till at the base is saturated with water, i.e. the
  till water amount is equal to `hydrology.tillwat_max`. This change should improve
  grounding line movement and make the basal yield stress modification turned on with
  `basal_yield_stress.slippery_grounding_lines` unnecessary.
- Fix the approximation of the driving stress at floating ice margins. (This fix was
  contributed by Ronja Reese and Torsten Albrecht.)
- Implement a mechanism for "optimizing" the till friction angle used by the Mohr-Coulomb
  yield stress model. The implementation is based on the code contributed by T. Albrecht.
- Add `atmosphere.elevation_change.precipitation.temp_lapse_rate` to the `-atmosphere
  ...,elevation_change` modifier. Now this parameter is used to compute the temperature
  change `dT` used to scale precipitation by a factor `exp(C * dT)` with `C =
  atmosphere.precip_exponential_factor_for_temperature`. We need a separate temperature
  lapse rate parameter to be able to use this modifier with atmosphere models that include
  an elevation-dependent near-surface air temperature parameterizations, e.g. `-atmosphere
  pik,elevation_change`.
- Improve the approximation of the grounding line flux (scalar and 2D diagnostics
  `grounding_line_flux`): the flux *through* the grounding line should be zero if its
  direction is parallel to the grounding line. Unfortunately this is impossible to achieve
  for an arbitrary grounding line shape if the grounding line is approximated by a mask on
  a uniform grid (as in PISM). This change improves the approximation for some
  combinations of grounding line shapes and grid resolutions. (This issue was reported by
  Ronja Reese.)

Changes from v1.1 to v1.2
=========================

Front retreat
^^^^^^^^^^^^^

- Implement the ISMIP6 front retreat parameterization. Reads a time-dependent ice extent
  mask (variable name: `land_ice_area_fraction_retreat`) from a file specified using the
  configuration parameter `geometry.front_retreat.prescribed.file`. This mechanism
  replaces the old "`ocean_kill`" calving code.
- Rename configuration parameters controlling front retreat because they are not
  calving-specific (`calving.front_retreat.use_cfl` to `geometry.front_retreat.use_cfl`
  and `calving.front_retreat.wrap_around` to `geometry.front_retreat.wrap_around`).
  Corresponding command-line options are renamed to `-front_retreat_cfl` and
  `-front_retreat_wrap_around`.
- Implement 3 frontal melt models: constant in time and space (`constant`), reading
  frontal melt from a file (`given`), and using time- and space-dependent thermal ocean
  forcing and modeled subglacial water flux in an implementation of the frontal melt
  parameterization in Rignot et al 2016 (`routing`).
- Now PISM combines retreat rates due to calving (eigen-calving and von Mises calving) and
  a frontal melt parameterizations before using these to update ice geometry. This
  simplifies and fixes the implementation of `geometry.front_retreat.use_cfl`.
- Add a configuration parameter `frontal_melt.include_floating_ice`: `true` means "apply
  computed frontal melt rates at *both* grounded and floating ice margins", `false` means
  "apply computed frontal melt rates at grounded margins only."

Subglacial hydrology
^^^^^^^^^^^^^^^^^^^^

- Add a new subglacial hydrology model `steady`. It adds an approximation of the
  subglacial water flux to the `null` model. This approximation uses the assumption that
  the subglacial water system instantaneously reaches the steady state after a change in
  the water input rate from the surface or a change in ice geometry. At high grid
  resolutions (~1km and higher) this is likely to be cheaper than using the `routing`
  model to obtain the flux needed by the frontal melt parameterization (above). Note that
  *this is a different model* and so when switching to it *re-tuning of the frontal melt
  parameterization will be necessary*.
- Variable `water_input_rate` containing water input rate for hydrology models (see
  `hydrology.surface_input.file`) uses units of "kg m-2 s-1" instead of "m s-1".
- Rename `hydrology.surface_input_file` to `hydrology.surface_input.file`. Also, add
  `hydrology.surface_input.period` and `hydrology.surface_input.reference_year` to support
  periodic water input rates.
- Add the ability to pass surface runoff modeled by PISM to a subglacial hydrology model
  (see `hydrology.surface_input_from_runoff`).
- Add a configuration parameter `hydrology.routing.include_floating_ice` to allow routing
  of subglacial water under floating ice. This may be appropriate when an outlet glacier
  has a small floating tongue. This also produces an extension of
  `subglacial_water_flux_mag` to floating areas, which is needed to model frontal melt
  using the `routing` model (see above).
- Add a configuration parameter `hydrology.add_water_input_to_till_storage`
  (default: yes). If "yes", surface water input is added to till storage which (if it
  overflows) then contributes to the amount of transportable subglacial water. If "no",
  surface water input directly contributes to the amount of transportable water, bypassing
  till storage. Basal melt rate always contributes to the amount of water stored in the
  till.

Surface, atmosphere, ocean forcing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Use linear interpolation in time for 2D time-dependent forcings that are interpreted as
  "snapshots" of a quantity. For example: ice surface temperature is interpreted as a
  snapshot while climatic mass balance is interpreted as a time-average over a specified
  interval. (Flux forcing fields such as the SMB are interpreted as piecewise-constant in
  time to simplify mass and flux accounting.) In particular, this change ensures that the
  2D sea level forcing results in a smoothly changing sea level.
- Fix units of the precipitation lapse rate (`(kg m-2/year)/km` instead of `(m/year)/km`).
- PISM uses configuration parameters to select surface, atmosphere, ocean, sea level, and
  frontal melt models. See `surface.models`, `atmosphere.models`, `ocean.models`,
  `sea_level.models`, and `frontal_melt.models`.
- Implement orographic precipitation following Smith and Barstad, *A linear theory of
  orographic precipitation*, 2004.
- Add regression tests for most of PISM's `surface`, `atmosphere`, `ocean`, and
  `sea_level` components.
- Rename `ocean.constant.melange_back_pressure_fraction` to
  `ocean.melange_back_pressure_fraction` and document its interaction with `-ocean
  ...,frac_MBP`.
- Rename `-atmosphere paleo_precip` to `precip_scaling`. Precipitation scaling using air
  temperature offsets is useful in other contexts, not just paleo-climate runs.
- Rename `-atmosphere lapse_rate` to `elevation_change`. This modifier includes
  adjustments that depend on the change in surface elevation but are *not* lapse rates.
- Rename `atmosphere.elevation_change.precipitation_lapse_rate` to
  `atmosphere.elevation_change.precipitation.lapse_rate` ("." instead of "_").
- Add `atmosphere.elevation_change.precipitation.method` (option `-precip_adjustment`).

  Set it to "scale" to use scaling with `exp(C * dT)`, where

  `C = atmosphere.precip_exponential_factor_for_temperature`

  and `dT` is computed using provided reference surface elevation and
  `atmosphere.elevation_change.temperature_lapse_rate`.

  Set it to "shift" to use `atmosphere.elevation_change.precipitation.lapse_rate` instead.
- Rename `-surface lapse_rate` to `-surface elevation_change`: now this modifier includes
  an adjustment that is not a lapse rate. Also, rename all related configuration
  parameters.
- Add a configuration parameter `surface.elevation_change.smb.method`: "scale" to use the
  exponential factor (see `surface.elevation_change.smb.exp_factor`); "shift" to use the
  lapse rate (see `surface.elevation_change.smb.lapse_rate`).
- Add `surface.elevation_change.smb.exp_factor`:

  `SMB = SMB_input * exp(C * dT),`

  where `C = surface.elevation_change.smb.exp_factor` and `dT` is the change in surface
  temperature computed using `surface.elevation_change.temperature_lapse_rate`.

Calving
^^^^^^^

- Add configuration parameters `calving.vonmises_calving.use_custom_flow_law`,
  `calving.vonmises_calving.flow_law`, and `calving.vonmises_calving.Glen_exponent`. This
  makes it possible to use a different flow law (and so a different Glen exponent) in the
  stress balance model and the tensile von Mises stress computation.

Energy balance
^^^^^^^^^^^^^^

- Switch to an unconditionally-stable method for the approximation of the heat equation in
  columns of the bedrock thermal layer (backward Euler time discretization instead of
  explicit time stepping).
- Add a configuration parameter `energy.bedrock_thermal.file`. Use this to specify a
  separate file containing the geothermal flux field (`bheatflx`). Leave it empty to read
  `bheatflx` from the input file (`-i` and `input.file`).

Bed deformation
^^^^^^^^^^^^^^^

- Fix the implementation of the elastic part of the Lingle-Clark bed deformation model.
  See `issue 424`_. To update this model we need to compute the discrete convolution of
  the current load with the load response matrix (the discrete equivalent of computing the
  convolution of the load with the Green's function).

  The load response matrix itself is approximated using quadrature and a one-dimensional
  interpolant for the tabulated Green's function.

  The old code used a "naive" implementation of the discrete convolution which was *both*
  slow and broken. The new implementation uses FFT to compute the discrete convolution,
  making it faster and (surprisingly) easier to implement.

  PISM now includes a regression test covering this. Unfortunately we don't have an exact
  solution to compare to, so the best we can do is this: a) compare PISM's FFT-based
  convolution to `scipy.signal.fftconvolve()` and b) compare PISM's load response matrix
  to an independent implementation using `scipy.integrate.dblquad()` (instead of
  `adapt_integrate()` by Steven G. Johnson).
- The configuration parameter `bed_deformation.lc.elastic_model` is set to "on" by
  default. This means that now `-bed_def lc` enables *both* the elastic and the viscous
  part of the Lingle-Clark model. In previous PISM versions `-bed_def lc` turned on the
  viscous part of the model and an extra command-line option (`-bed_def_lc_elastic_model`)
  was required to turn on the elastic part.
- Rename `bed_deformation.update_interval` to `bed_deformation.lc.update_interval` and fix
  its interpretation: before this change both bed deformation models (point-wise isostasy
  and the Lingle-Clark model) updated *not more often than* every
  `bed_deformation.update_interval` years. This lead to issues with stopped and re-started
  simulations (see `issue 422`_). Now the point-wise isostasy model is updated every time
  step (its computational cost is negligible) and the Lingle-Clark model is updated
  *exactly* every `bed_deformation.lc.update_interval` years, limiting PISM's time step
  length.

Regional modeling
^^^^^^^^^^^^^^^^^

- Improve handling of domain edges in regional simulations.
- Fix PISM's `-regional` runs: disable ice flow, surface mass balance, and basal mass
  balance effects on ice geometry in "no model" areas.

Stress balance
^^^^^^^^^^^^^^

- Improve PISM's handling of ice margins at locations where an icy cell is next to an
  ice-free cell that has surface elevation exceeding the surface elevation of ice (valley
  glaciers, fjords, nunataks). Here we have to use one-sided finite differences to avoid
  the influence of the surface elevation at ice-free locations on the estimate of the
  driving stress. In the SSAFD code we also prescribe drag along the margin (see
  `basal_resistance.beta_lateral_margin`) and a stress boundary condition assuming that at
  these locations ice is in hydrostatic equilibrium with the rock next to it. In other
  words, this boundary condition is the same as the calving front BC, but with a zero
  pressure difference.
- Rename command-line options `-ssa_rtol` to `-ssafd_picard_rtol` and `-ssa_maxi` to
  `-ssafd_picard_maxi` to make it clear that they control Picard iterations.

Basal strength
^^^^^^^^^^^^^^

- Add the ability to use space- and time-dependent `delta` (minimum effective pressure on
  till as a fraction of overburden pressure) in the Mohr-Coulomb basal yield stress
  parameterization.
- Yield stress models can be time-dependent.
- Implement "regional" versions of all yield stress models (both Mohr-Coulomb and
  constant). Previous versions did not support constant yield stress in regional model
  configurations.

Diagnostics
^^^^^^^^^^^

- Add configuration parameters `output.ISMIP6` (follow ISMIP6 conventions),
  `output.use_MKS` (save output variables in MKS units), `output.ISMIP6_extra_variables`
  (list of fields to report when `-extra_vars ismip6` is given), and
  `output.ISMIP6_ts_variables` (list of scalar time series to report when `-ts_vars
  ismip6` is given). When `output.ISMIP6` is set PISM saves spatially-variable diagnostics
  at the beginning of the run (if requested).
- Bug fix: ensure that land ice area fraction (diagnostic variable `sftgif`) never
  exceeds 1.
- Add the `hydraulic_potential` diagnostic to `routing` and `distributed` subglacial
  hydrology models.
- Fix an unreported bug in the computation of the `flux` diagnostic. This bug affected
  PISM's diagnostic variables `flux`, `velbar`, `velbar_mag`, and `vonmises_stress` (which
  uses `velbar`).

  It did *not* affect ice dynamics.
- Implement 2D and scalar grounding line flux diagnostics `grounding_line_flux`. See
  `issue #300`_.
- Rename `ocean_pressure_difference` to `ice_margin_pressure_difference`.

Input and output
^^^^^^^^^^^^^^^^

- Add a configuration flag `output.extra.stop_missing` (default: yes). Set this flag to
  "no" to make PISM warn about `-extra_vars` diagnostics that are requested but not
  available instead of stopping with an error message.
- Remove `output.variable_order`. Now PISM always uses `y,x,z` in output files.
- Allow using NCAR's ParallelIO library to write output files. If ParallelIO is compiled
  with parallel NetCDF-4 and PnetCDF libraries, this adds four new choices for
  `output.format`: `pio_pnetcdf` (parallel, using the CDF5 version of the NetCDF format),
  `pio_netcdf4p` (parallel, using the HDF5-based NetCDF format), `pio_netcdf4c` (serial,
  HDF5-based compressed NetCDF-4 format), `pio_netcdf` (serial).
- Fix `issue 327`_: now PISM uses mid-points of reporting intervals when saving to the
  `-extra_file`. This makes PISM's output files easier to process using CDO.

Other
^^^^^

- Provide better error messages when trying to read in a 2D field but the input file
  contains a 3D variable (or trying to read a 3D field but the input contains a 2D
  variable).
- Provide better error messages when trying to allocate more than 10000 records of a
  forcing field.
- PISM supports CMake 3.1 again (v1.1 required CMake 3.13 for no good reason).
- Use `pism/pism_config.hh` and `pism_config.cc` that are generated by CMake to include
  build information in the binary. This should take care of `issue 405`_.
- PISM uses the new (v5.x) PROJ API, so PROJ 5.0 or later is required to compute
  longitude-latitude grid coordinates and cell bounds. (Tested using PROJ v5.2.0 and
  v6.1.1.)
- Add contributing guidelines to the User's Manual.


Changes from v1.0 to v1.1
=========================

- PISM no longer attempts to use projection information to compute cell areas. This change
  was prompted by better mass accounting: it is now clear that using numerical methods
  designed for a uniform grid forces us to treat cells are equal in area. We may address
  this issue later but do not have the resources to work on this topic in the near future.
  Please use an equal-area projection for your simulations if distortions caused by
  working in a projected coordinate system are a concern.
- Add 5 more parameterizations of near-surface air temperature to `-atmosphere pik`.
- PISM stops with an error message if the name of a parameter in a `-config_override` file
  does not match any of the known PISM parameters.
- Fix `issue 375`_ (could not use `-config_override` to control the
  bed-elevation-dependent parameterization of the till friction angle).
- PISM stops with an error message if the diffusivity of the SIA flow exceeds a given
  threshold (see `stress_balance.sia.max_diffusivity`). Extremely high SIA diffusivities
  often mean that the setup is not "shallow enough"; in a situation like this it might
  make sense to re-evaluate model parameters before proceeding. (A short "smoothing" run
  might be helpful, too, if high diffusivities occur at the beginning of a simulation
  using ice thickness or bed topography not computed by PISM.)
- The SIA stress balance model limits computed diffusivity at
  `stress_balance.sia.max_diffusivity` if
  `stress_balance.sia.limit_diffusivity` is set. This makes it
  possible to speed up simulations in which high diffusivities at a
  few isolated grid points force PISM to take very short time steps.
  *This implies sacrificing accuracy at these grid points. Use with
  caution!*
- The SSAFD solver limits ice speed at a threshold specified by
  `stress_balance.ssa.fd.max_speed`. This may be useful when the computed sliding speed is
  abnormally high at a few isolated grid points, reducing the length of time steps PISM
  can take. Capping ice speed makes it possible to ignore troublesome locations and speed
  up some simulations. The default (500 km/year) is set high enough to deactivate this
  mechanism.
- Discard requested snapshot times that are outside of the modeled time interval. (This
  keeps PISM from overwriting a snapshot file written by one of the previous runs in a
  re-started simulation.)
- Add a new configuration parameter `stress_balance.sia.bed_smoother.theta_min` for the
  bed roughness parameterization in the SIA stress balance model.
- Added PICO, the *Potsdam Ice-shelf Cavity mOdel* (https://doi.org/10.5194/tc-2017-70).
  Use `-ocean pico` to enable and see the documentation of PISM's `ocean models`_ in the User's
  Manual for details.
- Added `-ocean ...,anomaly`, an ocean model *modifier* that reads spatially-variable
  sub-shelf mass flux anomalies from an input file.
- Exclude ice shelves from the ocean load provided to bed deformation models. See `issue
  363`_.
- Revert the change from v0.7 to v1.0 in the handling of energy conservation near ice
  margins. PISM v0.7 and earlier ignored horizontal advection and strain heating terms in
  the energy balance equation at grid points with neighbors below a given threshold ice
  thickness. PISM v1.0 eliminated this adjustment at ice margins. This version restores
  it, with the following additions.

  Set `energy.margin_ice_thickness_limit` to control
  the thickness limit used to trigger the special marginal treatment. Set parameters

  - `energy.margin_exclude_horizontal_advection`,
  - `energy.margin_exclude_vertical_advection`, and
  - `energy.margin_exclude_strain_heating`

  to control individual parts of this modification.

  The underlying issue is that the gradient of the ice thickness is discontinuous at
  grounded ice margins, and so its finite-difference approximation used by PISM is of poor
  quality. (The same applies to the gradient of the top surface elevation.) Errors in
  these approximations propagate to other quantities computed by PISM, notably the ice
  velocity and the strain heating. The poor quality of approximation of *these* quantities
  is the main reason for excluding them from the energy-balance computation.

  Preliminary tests show that excluding the strain heating term near ice margins is the
  most important modification.
- Fix `issue 400`_ (`viscous_bed_displacement` should not use the `coordinates`
  attribute).
- Support older (< 1.7) OpenMPI versions.
- Add a work-around needed to use old-ish NetCDF (4.0 - 4.1) with OpenMPI.
- Fix `issue 222`_ (`-part_grid` residual redistribution code used to lose mass in
  parallel runs).
- Add `geometry.part_grid.max_iterations` and increase it to 10.
- Fix a bug in `pismr -regional` (stored surface elevation was not initialized correctly)
- PDD model: add scalar diagnostics

  - `surface_accumulation_rate`,
  - `surface_melt_rate`,
  - `surface_runoff_rate`.

  See `issue 394`_. Also, rename `saccum`, `smelt`, `srunoff` to
  `surface_accumulation_flux`, `surface_melt_flux`, `surface_runoff_flux`
  respectively. Now PDD's climatic mass balance can be compared to the effective climatic
  mass balance: use `surface_accumulation_flux - surface_runoff_flux`.

  To save all these, use `-extra_vars` shortcuts `pdd_fluxes` and `pdd_rates`.
- PDD model: replace command-line options `-pdd_rand`, `-pdd_rand_repeatable` with one
  configuration parameter: `surface.pdd.method` (select from `expectation_integral`,
  `repeatable_random_process`, `random_process`).
- Fix `issue 74`_. (Now `basal_mass_flux_floating` is zero with the `float_kill`
  calving mechanism, i.e. when `ice_area_glacierized_floating` is zero.)
- Refactor hydrology models, adding proper mass accounting.
- Implement 2D diagnostics quantities needed for mass conservation accounting in hydrology
  models:

  - `tendency_of_subglacial_water_mass`,
  - `tendency_of_subglacial_water_mass_due_to_input`,
  - `tendency_of_subglacial_water_mass_due_to_flow`,
  - `tendency_of_subglacial_water_mass_due_to_conservation_error`,
  - `tendency_of_subglacial_water_mass_at_grounded_margins`,
  - `tendency_of_subglacial_water_mass_at_grounding_line`, and
  - `tendency_of_subglacial_water_mass_at_domain_boundary`.

  Use the shortcut `hydrology_fluxes` to save all these in an "extra file."
- Add `hydrology.surface_input_file`: `IceModel` can read in time-dependent 2D water
  input rates for subglacial hydrology models.
- Implement a proper generalization to 2D of the 1D parameterization of the grounding line
  position. (This code interprets ice thickness, bed elevation, and sea level as
  piecewise-linear functions on a specially-designed triangular mesh refining the regular
  grid used by PISM.)
- Support 2D (spatially-variable) sea level elevation everywhere in PISM, including 2D sea
  level forcing. (Use `-sea_level constant,delta_sl_2d` and search for
  `ocean.delta_sl_2d.file` and related configuration parameters.)
- Split sea level forcing from the ocean model so that the sea level is available when
  sub-shelf melt parameterizations are initialized. Use `-sea_level constant,delta_sl`
  instead of `-ocean constant,delta_SL`.
- Decouple calving law parameterization from ocean models and the stress balance code.
- Add regression tests for all ocean models.
- Fix `issue 402`_: ensure reproducibility of `-bed_def lc` results.
- Clean up PISM's ocean, surface, and atmosphere model code, making it easier to test and
  debug.
- Make it easier to use scalar and 2D time-dependent forcing fields.
- Add configuration parameters `input.file` and `input.bootstrap`, corresponding to
  command-line options `-i` and `-bootstrap`.
- Add notes documenting the implementation of the calving front boundary condition to the
  manual.
- Make it easier to "balance the books":

  #. rename scalar diagnostics so that they match 2D diagnostics and
  #. report fluxes in `Gt/year` instead of `kg/year`.
- Update the Debian/Ubuntu section of the installation manual.
- Move the documentation of the BOMBPROOF numerical scheme for energy conservation from
  the source code browser into the manual.
- Add an experimental implementation of a parameterization of cryo-hydrologic warming
  based on *Cryo-hydrologic warming: A potential mechanism for rapid thermal response of
  ice sheets* by Phillips et al, 2010.)

Changes from v0.7 to v1.0
=========================

This document lists notable changes from PISM v0.7 to v1.0.

Summary
-------

- New mass transport code makes it easier to "balance the books".
- PISM's grids are no longer transposed ( ``(y,x)`` versus ``(x,y)`` ).
- Adds an optimized implementation of the GPBLD flow law for the Glen n=3 case.
- Adds von Mises calving (see Morlighem et al, *Modeling of Store Gletscher's calving
  dynamics, West Greenland, in response to ocean thermal forcing*, 2016)
- Adds more diagnostic quantities (127 spatially-variable fields and 38 scalar variables
  in total)
- Better code, `better documentation`_, more regression and verification tests.

Please run ``git log v0.7..v1.0`` for the full list.

See files in the ``doc/`` sub-directory for changes from v0.6 to v0.7, etc.

Installation
------------

- Remove ``Pism_BUILD_TYPE`` and use ``CMAKE_BUILD_TYPE`` instead.

Prerequisites
^^^^^^^^^^^^^

- Require CMake 3.1 and compilers supporting C++11.

- Require PETSc built with ``PetscScalar`` as ``double``. Stop if ``PetscScalar`` is
  ``complex``. See `issue 237`_.

- Drop Subversion support. Please use Git to download PISM source code.

- PETSc < 3.5 is not supported; use PETSc 3.5 and newer (PETSc 3.6.0 is not supported due
  to a bug).

Library and directory structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Install PISM headers in ``include/pism``, skipping 3rd party headers and empty
  directories (see `issue 166`_.)

- Link all of PISM into one single library.

- Install all Python scripts in ``util/``. Fixes `issue 346`_.

- Fix the directory structure created by ``make install``.

Other
^^^^^

- Remove all ``simpleXXX`` executables. See `issue 343`_. Use Python wrappers to access exact
  solutions used in PISM's verification tests.

- Remove ``pismo`` (use ``pismr -regional``).

Documentation
-------------

- Migrate documentation to Sphinx_.

- New PISM support e-mail address: uaf-pism@alaska.edu instead of help@pism-docs.org.

Computational domain and grid
-----------------------------

- Add options ``-x_range``, ``-y_range``, which specify domain extent in the `x` and `y`
  direction during bootstrapping. These can be used to extract a subset of a grid for a
  regional run.

- De-couple grid periodicity from grid registration and add the ``grid.registration``
  parameter. This changes the interpretation of ``-Lx`` and ``-Ly`` during bootstrapping.
  See `issue 347`_.

- Support EPSG:26710, EPSG:3413, and EPSG:3031. When an input file contains the global
  attribute ``proj4`` containing the string "``+init=epsg:XXXX``" where ``XXXX`` is one of
  these codes PISM will create a CF-conforming ``mapping`` variable with projection
  parameters corresponding to the selected mapping. See `issue 350`_.

- Write PROJ.4 parameters to ``mapping:proj4_params`` (for CDO).

Ice rheology
------------

- Add ``gpbld3``, the ``n==3`` optimized flow law.

  This is an optimized (vectorized_) implementation of the
  Glen-Paterson-Budd-Lliboutry-Duval flow law with the fixed Glen exponent of 3.

  On modern (2011 and on) CPUs this flow law implementation is almost 4 times faster than
  the default one. This significantly reduces the cost of high-resolution runs.

  The implementation uses ``exp()`` from VDT_, a vectorized math library developed at CERN.
  To reduce the number of external dependencies a copy of VDT (v0.3.6) is included in
  PISM's source tree.

Stress balance
--------------

- SSAFD KSP solver: use the initial residual norm.

  This prevents the SSAFD solver from failing when the solver has no work to do.

- Make the SSAFD solver a little more robust by replacing zero diagonal matrix entries
  with large beta, effectively "disabling" sliding at these locations. See `issue 349`_.

- Remove ``SIA_Sliding``, EISMINT II tests G and H, verification test E.

- Add ``stress_balance.vertical_velocity_approximation``. I.e. (optionally) use
  first-order upwinding to compute u_x and v_y in the vertical velocity computation.

- Add enhancement factors for interglacial periods (See Ralf Greve, *Application of a
  polythermal three-dimensional ice sheet model to the Greenland ice sheet: Response
  to steady-state and transient climate scenarios*, 1997.)

  Use the following configuration parameters to control this:

  - ``stress_balance.sia.enhancement_factor_interglacial``
  - ``stress_balance.ssa.enhancement_factor_interglacial``
  - ``time.eemian_start``
  - ``time.eemian_end``
  - ``time.holocene_start``

Geometry and mass transport
---------------------------

- Completely redesign and re-implement the mass transport code. The new code is
  well-isolated and extensible, designed to make "balancing the books" easier, and can be
  tested in isolation. See also `issue 201`_.

- Add the class ``Geometry`` that can be used to provide geometry information to PISM's
  sub-models. This improves interfaces of PISM's sub-models, reducing undesirable "tight"
  coupling.

- Option ``-part_grid`` implies ``-part_redist``.

Calving
-------

- Generalize eigen-calving code and add von Mises calving.

- Implement calving front retreat due to frontal melting.

- Rename ``-cfl_eigen_calving`` to ``-calving_cfl``.

- Make it possible to disable ``float_kill`` near grounding lines. See
  ``-float_kill_calve_near_grounding_line``.

- Add option ``-float_kill_margin_only``. See `issue 340`_.

- Allow using spatially-variable calving at thickness thresholds.

- Add ``-calving_wrap_around`` for synthetic geometry setups.

Energy conservation
-------------------

- ``BedThermalUnit`` ensures that computed bedrock temperatures exceed
  zero kelvin. See `issue 313`_.

- PISM no longer ignores horizontal enthalpy advection and strain
  heating near ice margins. See `issue 292`_.

- Following a re-interpretation of Aschwanden et al, *An enthalpy formulation for glaciers
  and ice sheets*, 2012 we require that dH/dp=0.

  Assuming that specific heat capacities of ice and water do not depend on temperature,
  this gives

  ``L(p) = (T_m(p) - T_m(p_air)) (c_w - c_i) + L_0``, where

  .. csv-table::

     ``T_m``   , melting temperature
     ``c_w``   , specific heat capacity of water
     ``c_i``   , specific heat capacity of ice
     ``L_0``   , latent heat of fusion at air pressure
     ``p_air`` , air (atmospheric) pressure

  Note that this form of the latent heat of fusion ``L(p)`` also follows from Kirchhoff's
  law of thermochemistry. See ``EnthalpyConverter::L(T_pm)`` for details. See `issue
  334`_.

- To allow for better code optimization, ``EnthalpyConverter`` no longer uses virtual
  methods. ``ColdEnthalpyConverter`` used in temperature-based verification tests sets ice
  melting temperature to 1e6 kelvin to ensure that all ice is considered "cold."
  ``varcEnthalpyConverter``, which implemented linear-in-temperature specific heat
  capacity of ice, is removed.

- Code solving the enthalpy equation within an ice column supports both Dirichlet and
  Neumann boundary conditions at the top surface.

  Only the Dirichlet condition is used in modeling runs; Neumann B.C. code is there to
  simplify testing.

- Documented the discretization of the enthalpy column system. Added simple verification
  tests for the enthalpy solver within an ice column (pure advection and pure diffusion
  with different boundary conditions).

- To simplify model initialization and testing energy balance models are isolated. The
  rest of PISM uses the interface class ``EnergyModel``. The old "cold mode"
  temperature-based energy balance model is in ``TemperatureModel``. The enthalpy-based
  model is in ``EnthalpyModel``.

Input and output
----------------

- Remove the HDF5-based parallel I/O code.

- Remove ``-o_format quilt`` and ``pismmerge``.

- Implement reading string attributes from NetCDF-4 files.

- Add detailed I/O (writing) reporting with ``-verbose 3``.

- Add ``pism::StringLogger``, a logger that prints to a string.

- Add an option ``-profile`` to write detailed profiling information.

- Add ice thickness thresholds for reporting and stress balance.

  This makes it easier to track changes corresponding to "glacierized" areas while
  excluding the seasonal cover.

  See ``output.ice_free_thickness_standard`` and
  ``stress_balance.ice_free_thickness_standard``.

- Write run statistics to extra and time-series files. (See `issue 324`_, `issue 330`_.)

- New option: ``-save_force_output_times``.

- Avoid re-writing metadata that does not change during the run.

Diagnostics
^^^^^^^^^^^

- Add numerous new diagnostic quantities, including sets of diagnostics needed to "balance
  the books" when accounting for mass changes (conservation).

- Add scalar diagnostics using the new (higher) thickness threshold used to determine if a
  cell ice "ice-free". These diagnostics have the "``_glacierized``" suffix and can be
  interpreted as tracking changes in glacierized areas (ignoring the seasonal cover).

- Rates of change reported by PISM are *mean* rates of change over reporting intervals
  computed using finite differences.

- Better feedback on missing (or renamed) diagnostics. If a requested diagnostic is not
  available PISM will stop with an error message listing available diagnostics.

Bed deformation
---------------

- Add a new command-line option: ``-uplift_file``. Use it to specify the name of a file
  containing the variable ``dbdt`` to use when initializing the Lingle-Clark bed
  deformation model. See `issue 390`_.

- Add ``-topg_delta_file topg_delta.nc.``

  With this option PISM tries to read "topg_delta" from a specified file and sets bed
  topography at the beginning of a run to

  .. code::

     bed_elevation = topg + topg_delta.

  Here ``topg`` is read from an input file (``-i``), ``topg_delta`` -- from
  ``topg_delta.nc``.

- Lingle-Clark bed deformation model: save the viscous bed displacement on the extended
  grid so that stopping and re-starting the model does not affect results. This also makes
  it possible to refine computational grids in runs using the model. See `issue 370`_.

- Bed deformation models can be used and tested in isolation (see `issue 181`_).

Subglacial hydrology
--------------------

- Re-implement lateral till water diffusion as in Bueler and Brown, 2009.

Climate forcing
---------------

- Apply lapse rate corrections throughout the domain.

  Previously it was used in icy areas only.

- Remove old PDD code.

- ``-atmosphere``: use "``kg m-2 second-1``" precipitation units.

- Add ``ocean_frac_SMB``, a modifier scaling shelf-base mass flux

- Atmosphere and ocean modifiers save "effective" fields.

- Add an option and config. parameter ``surface.force_to_thickness.start_time`` to allow
  delaying the nudging effect.

Bug fixes
---------

(This is an incomplete list.)

- Fix `issue 328`_ (diagnostic computation of ``wvelsurf``).

- Fix a bug in ``pism::ocean::Constant`` (``-shelf_base_melt_rate`` was ignored).

- Fix `issue 351`_ (duplicate history in -extra_files).

- Fix a bug in the code implementing ``-save_file`` with ``-save_split`` (see `issue 325`_).

- Fix `issue 323`_ (fix EISMINT II settings so v0.7 conforms).

- Fix `issue 321`_: Sea level affects margin stress B.C. in the "dry simulation" mode.

- Fix interpolation weights and add a test. See `issue 326`_.

Miscellaneous
-------------

- Undo the "fundamental transpose": now PISM uses the (y,x) order in files and memory.

  This simplifies pre-processing of input files and post-processing and analysis of
  modeling results.

- Allow extrapolation during regridding to simplify restarting in runs where ice thickness
  exceeded the height of the computational domain *and* to extend the domain in
  continental ice sheet simulations. See `issue 302`_.

- Save the model state if the ice thickness exceeds the height of the computational
  domain.

- The age model was moved to ``AgeModel``.

- Add the ability to add "hooks" to ``RuntimeError``.

  Added to allow custom actions (such as printing a traceback) when an error is detected.

- Improve PISM's version information

  - Add committer's name and date to the version string.
  - ``pismr -version`` prints versions of

    - PISM
    - PETSc (including configuration options)
    - MPI
    - NetCDF
    - FFTW
    - GSL
    - PROJ.4
    - SWIG (if Python bindings are enabled)

- Add support for coverage testing using ``lcov``.

  Set ``Pism_CODE_COVERAGE`` to enable, use ``make coverage_report`` to generate a report and
  and ``make coverage_reset`` to reset coverage data.

- Add ``.clang-format`` to the top level directory

  ``clang-format`` makes it much easier to use consistent code formatting throughout. To
  re-format a file, commit it to the repository, then run

  ``clang-format -i filename.cc``

  (Here ``-i`` tells clang-format to edit files "in place." Note that editing in place
  is safe because you added it to the repository.)

- Re-organize configuration parameters: all parameters have new names that reflect their
  places within the model hierarchy.

- Improve processing of boolean command-line options

  .. code::

     -foo yes
     -foo on
     -foo true
     -foo True
     -foo (no argument)

  set the boolean flag to "true."

  .. code::

     -foo no
     -foo false
     -foo False
     -no_foo (for backward compatibility)

  set the flag to "false."

- Add numerous regression tests.

.. _Sphinx: http://pism.io/docs
.. _better documentation: Sphinx_
.. _vectorized: https://en.wikipedia.org/wiki/Automatic_vectorization
.. _VDT: https://github.com/dpiparo/vdt

.. _issue 74:  https://github.com/pism/pism/issues/74
.. _issue 166: https://github.com/pism/pism/issues/166
.. _issue 181: https://github.com/pism/pism/issues/181
.. _issue 201: https://github.com/pism/pism/issues/201
.. _issue 222: https://github.com/pism/pism/issues/222
.. _issue 237: https://github.com/pism/pism/issues/237
.. _issue 292: https://github.com/pism/pism/issues/292
.. _issue 300: https://github.com/pism/pism/issues/300
.. _issue 302: https://github.com/pism/pism/issues/302
.. _issue 313: https://github.com/pism/pism/issues/313
.. _issue 321: https://github.com/pism/pism/issues/321
.. _issue 323: https://github.com/pism/pism/issues/323
.. _issue 324: https://github.com/pism/pism/issues/324
.. _issue 325: https://github.com/pism/pism/issues/325
.. _issue 326: https://github.com/pism/pism/issues/326
.. _issue 327: https://github.com/pism/pism/issues/327
.. _issue 328: https://github.com/pism/pism/issues/328
.. _issue 330: https://github.com/pism/pism/issues/330
.. _issue 334: https://github.com/pism/pism/issues/334
.. _issue 340: https://github.com/pism/pism/issues/340
.. _issue 343: https://github.com/pism/pism/issues/343
.. _issue 346: https://github.com/pism/pism/issues/346
.. _issue 347: https://github.com/pism/pism/issues/347
.. _issue 349: https://github.com/pism/pism/issues/349
.. _issue 350: https://github.com/pism/pism/issues/350
.. _issue 351: https://github.com/pism/pism/issues/351
.. _issue 370: https://github.com/pism/pism/issues/370
.. _issue 390: https://github.com/pism/pism/issues/390
.. _issue 394: https://github.com/pism/pism/issues/394
.. _issue 400: https://github.com/pism/pism/issues/400
.. _issue 402: https://github.com/pism/pism/issues/402
.. _issue 363: https://github.com/pism/pism/issues/363
.. _issue 405: https://github.com/pism/pism/issues/405
.. _issue 422: https://github.com/pism/pism/issues/422
.. _issue 424: https://github.com/pism/pism/issues/424
.. _issue 407: https://github.com/pism/pism/issues/407
.. _ocean models: http://www.pism.io/docs/climate_forcing/ocean.html
..
   Local Variables:
   fill-column: 90
   End:
