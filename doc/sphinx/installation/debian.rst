.. _sec-install-debian:

Installing prerequisites from Debian packages
---------------------------------------------

You should be able to use your package manager to get the prerequisites for PISM. Install
the following packages using ``apt-get`` or ``synaptic`` or similar. All of these are
recommended as they satisfy requirements for building or running PISM.

.. csv-table:: Debian packages
   :header: Package name, Comments

   ``cmake``,            required to configure PISM
   ``libfftw3-dev``,     required by PISM
   ``g++``,              required to build PISM
   ``libgsl-dev``,       required by PISM
   ``netcdf-bin``,       required: ``ncgen`` is used during the build process
   ``libnetcdf-dev``,    required by PISM
   ``libudunits2-dev``,  required by PISM
   ``cdo``,              used in some pre-processing scripts
   ``cmake-curses-gui``, a text-based easy interface for CMake
   ``git``,              used to get PISM source code
   ``nco``,              used in many pre-processing scripts
   ``ncview``,           view fields in NetCDF files
   ``libproj-dev``,      used to compute longitude and latitude coordinates of grid points
   ``python-dev``,       (helps with scripts…perhaps not essential)
   ``python-pyproj``,    used in some pre-processing scripts
   ``python-netcdf4``,   used in most post-processing scripts
   ``libx11-dev``,       X windows is useful to get graphics through PETSc
   ``libblas-dev``,      BLAS is required by PETSc
   ``liblapack-dev``,    LAPACK is required by PETSc
   ``mpich``,            MPI is required to run PISM in parallel
   ``libmpich-dev``,     MPI is required to run PISM in parallel

You may be able to install these by running

.. literalinclude:: code/install_libraries.sh
   :language: bash
   :lines: 3-

.. only:: html

   Click :download:`here <code/install_libraries.sh>` to download this file.

(You may need to change this command to match your package system.)

Once done, see :ref:`sec-install-petsc` to install PETSc from source and then
:ref:`sec-install-pism` for building PISM itself.
