.. include:: ../prologue.rst

.. _sec-install-prerequisites:

Required libraries
==================

This table lists required dependencies for PISM alphabetically.

.. csv-table::
   :header: Required Library, Comment

   FFTW_,        version >= 3.1
   GSL_,         version >= 1.15
   MPI_,         any recent version
   NetCDF_ [#]_, version >= 4.1
   PETSc_ [#]_,  version >= 3.5
   UDUNITS_,     any recent version

Before installing these "by hand", check sections :ref:`sec-install-debian` and
:ref:`sec-install-macos` for specific how-to.

In particular, if multiple MPI implementations (e.g. MPICH and Open-MPI) are installed
then PETSc can under some situations "get confused" and throw MPI-related errors. Even
package systems have been known to allow this confusion.

Optional libraries listed below are needed for certain PISM features, namely cell-area
correction and parallel I/O. These libraries are recommended, but not strictly required:

.. csv-table::
   :header: Recommended Library, Comment

   PROJ.4_,  Used to compute cell areas and cell bounds
   PnetCDF_, Can be used for parallel I/O

Python_ is needed both in the PETSc installation process and in scripts related to PISM
pre- and post-processing, while Git_ is usually needed to download the PISM code. Both
should be included in any Linux/Unix distribution.

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
