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
   * - ``util/fill_missing.py``
     - Uses an approximation to Laplace's equation :math:`\nabla^2 u = 0` to smoothly
       replace missing values in a two-dimensional NetCDF variable. The "hole" is filled
       with an average of the boundary non-missing values. Depends on ``netcdf4-python``
       and ``scipy`` Python packages.
   * - ``util/flowline.py``
     - See section :ref:`sec-flowline-modeling`.
   * - ``util/flowlineslab.py``
     - See section :ref:`sec-flowline-modeling`.
   * - ``util/check_stationarity.py``
     - Evaluate stationarity of a variable in a PISM ``-ts_file`` output.
   * - ``util/nc2cdo.py``
     - Makes a netCDF file ready for Climate Data Operators (CDO).
   * - ``util/nc2mat.py``
     - Reads specified variables from a NetCDF file and writes them to an output file in
       the MATLAB binary data file format ``.mat``, supported by MATLAB version 5 and
       later. Depends on ``netcdf4-python`` and ``scipy`` Python packages.
   * - ``util/nccmp.py``
     - A script comparing variables in a given pair of NetCDF files; used by PISM software
       tests.
   * - ``util/pism_config_editor.py``
     - Makes modifying or creating PISM configuration files easier.
   * - ``util/pism_matlab.m``
     - An example MATLAB script showing how to create a simple NetCDF file PISM can
       bootstrap from.
   * - ``util/PISMNC.py``
     - Used by many Python example scripts to generate a PISM-compatible file with the
       right dimensions and time-axis.
