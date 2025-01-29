This directory contains

- A script that generates time-dependent orbital parameters that can be used for
  paleo-simulations with dEBM-simple. It uses trigonometric expansions from

  A. L. Berger, "Long-Term Variations of Daily Insolation and Quaternary Climatic
  Changes," Journal of the Atmospheric Sciences, vol. 35, no. 12, pp. 2362–2367, Dec.
  1978.

  The implementation itself is heavily based on the code from GISS ModelE version 2.1.2
  (see https://simplex.giss.nasa.gov/snapshots/modelE2.1.2.tgz ).

  Note: ModelE documentation states "our code is freely available for download and use",
  but unfortunately it does not include a license.

  Run `orbital_parameters.py --help` for details.

  Note: PISM requires time bounds for all time-dependent inputs; see PISM's manual or the
  `Makefile` for a way to use `ncap2` from NCO to add them.

- A script `run.sh` and a `Makefile` that uses a SeaRISE-Greenland setup to run
  dEBM-simple

  - in the "present day" configuration
  - in a "paleo" configuration with time-dependent orbital parameters
  - in a "present day" configuration with provided time-dependent surface albedo

  and save a number of diagnostics to NetCDF files.

  Please see `run.sh` for command-line options needed by these simulations.
