.. default-role:: literal

Changes from v1.1.3 to v1.1.4
=============================

- PISM can be built with PROJ v6. (We define `ACCEPT_USE_OF_DEPRECATED_PROJ_API_H`. This
  workaround will break once PROJ drops the deprecated API completely. See `issue 409`_.)
- Fix a minor bug in the `routing` hydrology model (improper indexing in the code
  computing hydraulic conductivity).

Changes from v1.1.2 to v1.1.3
=============================

- Minor fixes of PISM's documentation.
- Fix an unreported bug in the computation of the `flux` diagnostic. This bug affected
  PISM's diagnostic variables `flux`, `velbar`, `velbar_mag`, and `vonmises_stress` (which
  uses `velbar`).

  It did *not* affect ice dynamics.

Changes from v1.1.1 to v1.1.2
=============================

- Fix an unreported bug in `-surface ...,forcing`: PISM was ignoring the time step
  restriction associated with this mechanism; large `surface.force_to_thickness.alpha`
  could lead to uncontrolled growth of ice thickness.
- Update the `-atmosphere pik` temperature parameterization for compatibility with paleo
  simulations by Albrecht et al.
- Switch to an unconditionally-stable method for the approximation of the heat equation in
  columns of the bedrock thermal layer (backward Euler time discretization instead of
  explicit time stepping).

Changes from v1.1 to v1.1.1
===========================

- PISM supports CMake 3.1 again (v1.1 required CMake 3.13 for no good reason).
- Fix PISM's `-regional` runs: disable ice flow, surface mass balance, and basal mass
  balance effects on ice geometry in "no model" areas.

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
  zero Kelvin. See `issue 313`_.

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
  melting temperature to 1e6 Kelvin to ensure that all ice is considered "cold."
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

.. _Sphinx: http://pism-docs.org/sphinx/
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
.. _issue 302: https://github.com/pism/pism/issues/302
.. _issue 313: https://github.com/pism/pism/issues/313
.. _issue 321: https://github.com/pism/pism/issues/321
.. _issue 323: https://github.com/pism/pism/issues/323
.. _issue 324: https://github.com/pism/pism/issues/324
.. _issue 325: https://github.com/pism/pism/issues/325
.. _issue 326: https://github.com/pism/pism/issues/326
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
.. _issue 409: https://github.com/pism/pism/issues/409
.. _ocean models: http://pism-docs.org/sphinx/climate_forcing/ocean.html
..
   Local Variables:
   fill-column: 90
   End:
