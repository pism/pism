.. include:: ../../global.txt

.. _sec-scripts:

Utility and test scripts
------------------------

In the ``test/`` and ``util/`` subdirectories of the PISM directory the user will find
some python scripts and one Matlab script, listed in :numref:`tab-scripts-overview`. The
python scripts are all documented at the *Packages* tab on the `PISM Source Code
Browser <pism-browser_>`_. The Python scripts all take option ``--help``.

.. list-table:: Some scripts which help in using PISM
   :name: tab-scripts-overview
   :header-rows: 1
   :widths: 1,1

   * - Script
     - Function
   * - ``test/vfnow.py``
     - Organizes the process of verifying PISM. Specifies standard refinement paths for
       each of the tests (section :ref:`sec-verif`).
   * - ``test/vnreport.py``
     - Automates the creation of convergence graphs like figures :numref:`fig-thickerrsB`
       -- :numref:`fig-velerrsI`.
   * - ``util/pism_fill_missing``
     - Uses an approximation to Laplace's equation :math:`\nabla^2 u = 0` to smoothly
       replace missing values in a two-dimensional NetCDF variable. The "hole" is filled
       with an average of the boundary non-missing values. Depends on ``netcdf4-python``
       and ``scipy`` Python packages.
   * - ``util/pism_flowline``
     - See section :ref:`sec-flowline-modeling`.
   * - ``util/pism_check_stationarity``
     - Evaluate stationarity of a variable in a PISM ``-ts_file`` output.
   * - ``util/pism_nc2cdo``
     - Makes a netCDF file ready for Climate Data Operators (CDO).
   * - ``util/pism_nccmp``
     - A script comparing variables in a given pair of NetCDF files; used by PISM software
       tests.
   * - ``examples/preprocessing/flowlineslab.py``
     - See section :ref:`sec-flowline-modeling`.
   * - ``examples/preprocessing/PISMNC.py``
     - Used by several Python example scripts to generate a PISM-compatible file with the
       right dimensions and time-axis.
