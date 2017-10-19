.. include:: ../global.txt

.. _sec-install-sources:

Installing prerequisites from sources
-------------------------------------

#. You will need Python_ and Git_ installed. To use the (recommended) graphical output of
   PISM you will need an `X Window server <X_>`_.

#. Generally the "header files" for its prerequisite libraries are required for building
   PISM. (This means that the "developer’s versions" of the libraries are needed if the
   libraries are downloaded from package repositories like Debian's; see
   :ref:`sec-install-prerequisites`.)

#. PISM uses `NetCDF <NetCDF_>`_ as an input and output file format. If it is not already
   present, install it using the instructions at the web-page or using a package
   management system.

#. PISM uses the `GNU Scientific Library <GSL_>`_ for certain numerical calculations and
   special functions. If it is not already present, install it using the instructions at
   the web-page or using a package management system.

#. PISM uses the `FFTW library <FFTW_>`_ for the deformation of the solid earth (bed)
   under ice loads. Install FFTW version 3.1 or later, or check that it is installed
   already.

#. You will need a version of `MPI <MPI_>`_. Your system may have an existing MPI
   installation, in which case it should probably be used when building PETSc. The goal is
   to have the PETSc installation use the same version of MPI which is called by the
   ``mpiexec`` or ``mpirun`` executable.

   If you had to install an MPI library "by hand" you will want to add the MPI ``bin``
   directory to your path so that you can run parallel programs using the ``mpiexec`` or
   ``mpirun`` command. For example, you can add it with the statement

   .. code-block:: bash

      export PATH=/home/user/mympi/bin:$PATH

   (for Bash shell).

   Such a statement can, of course, appear in your ``.bashrc`` (or
   ``.profile``) file so that there is no need to retype it each time
   you use MPI.

#. PISM uses UDUNITS_ to convert units of physical quantities read from input files and
   written to output files. Follow instructions on its website to install.
