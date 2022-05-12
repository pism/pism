.. include:: ../global.txt

.. _sec-install-prerequisites:

Required tools and libraries
============================

This table lists required dependencies for PISM alphabetically.

.. csv-table::
   :header: Required Library, Comment

   CMake_,       version `\ge` 3.3
   FFTW_,        version `\ge` 3.1
   GSL_,         version `\ge` 1.15
   MPI_,         any recent version
   NetCDF_ [#]_, version `\ge` 4.4
   PETSc_ [#]_,  version `\ge` |petsc-min-version|
   UDUNITS_,     any recent version

Before installing these "by hand", check sections :ref:`sec-install-debian` and
:ref:`sec-install-macos` for specific how-to.

In particular, if multiple MPI implementations (e.g. MPICH and Open MPI) are installed
then PETSc can under some situations "get confused" and throw MPI-related errors. Even
package systems have been known to allow this confusion.

.. note::

   We recommend un-installing all MPI libraries except one. In most cases there is no
   reason to have both Open MPI and MPICH, for example.

Optional libraries listed below are needed for certain PISM features, namely computing
longitude, latitude coordinates of grid points and parallel I/O. These libraries are
recommended, but not strictly required:

.. csv-table::
   :header: Optional Library, Comment

   PROJ_,  version `\ge` 6.0 (used to compute longitude-latitude grid coordinates and cell bounds)
   PnetCDF_, Can be used for faster parallel I/O
   ParallelIO_, Can be used for faster parallel I/O

Python_ is needed for the PETSc installation process; a number of PISM's pre- and
post-processing scripts also use Python (either 2.7 or 3.x), while Git_ is usually needed
to download the PISM code.

PISM's Python bindings support Python 2.7 or 3.3 and later [#]_.

The following Python packages are needed to do all the examples in the :ref:`User’s Manual
<sec-users-manual>` (which run Python scripts):

.. csv-table:: Python packages
   :header: Library, Comment

   NumPy_,          used in *most* scripts
   matplotlib_,     used in some scripts
   netcdf4-python_, used in *most* scripts

.. rubric:: Footnotes

.. [#] Note that PISM uses ``ncgen`` (provided by NetCDF) on the system where PISM is
       *compiled*.
.. [#] "PETSc" is pronounced "pet-see".
.. [#] PISM's Python bindings are tested using Python 3.6.
