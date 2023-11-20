Changes from 0.4 (June 2011) to 0.5
===================================

-  Switch to PETSc 3.2
-  Move to github.com
-  Add regional modeling tools: pismo, Python tools
-  Improve the IceFlowLaw class by moving physical constants out of it
-  Implement ``-o_format [netcdf4_parallel, pnetcdf]``
-  Implement ``-o_order [xyz, yxz, zyx]``
-  Implement ``-e_age_coupling``
-  Improved the implementation of PISM's basal yield stress model
-  Separated flow laws used in SIA and SSA code. New command-line
   options: ``-sia_e``, ``-ssa_e``. New config. parameters:
   ``sia_enhancement_factor``, ``ssa_enhancement_factor``
-  Fixed ``-topg_to_phi``: now it takes 4 numbers (no special value in
   the ocean)
-  ``-pdd_annualized`` implemented, and recommended for degree-day
   scheme usage
-  Flush time-series when ``-extra_files`` are written
-  Improved documentation (mostly climate forcing code)
-  Improved climate forcing code and its interface (in particular:
   consistent command-line options)
-  Fixed numerous bugs
-  Improved regridding code
-  Implemented temperature-dependent conductivity and specific heat
   capacity of ice
-  Model time is in seconds
-  PISM uses time bounds, both in reporting code and code using forcing
   data
-  Reported rates of change are computed as average rates over reporting
   intervals, by differencing cumulative quantities
-  Reporting: daily, monthly, yearly
-  Added support for the Gregorian calendar
-  PISM uses Proj.4 to compute latitudes, longitudes, and cell areas
-  Implemented ``-tauc_to_phi``
-  Updated and improved modeling examples
-  Removed the pgrn executable; all whole-Greenland examples use pismr
-  Removed the EISMINT-Ross example and the pross executable
-  Inverse modeling tools (not part of the release)
