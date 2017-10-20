.. include:: ../../global.txt

.. _sec-nctoolsintro:

Handling NetCDF files
---------------------

PISM takes one or more NetCDF files as input, performs some computation, and then produces
one or more NetCDF files as output. However, other tools are usually needed to help to
extract meaning from NetCDF files, and yet more NetCDF tools help with creating PISM input
files or post-processing PISM output files.

Here we list a number of NetCDF tools that can be useful in preparing input data for use
with PISM and post-processing results; see :numref:`tab-NetCDFview`.

.. list-table:: A selection of tools for viewing and modifying NetCDF files.
   :name: tab-NetCDFview
   :header-rows: 1

   * - Tool
     - Function

   * - ``ncdump`` (part of NetCDF_)
     - dump binary NetCDF as ``.cdl`` (text) file

   * - ``ncgen`` (part of NetCDF_)
     - convert ``.cdl`` file to binary NetCDF

   * - ncview_
     - quick graphical view

   * - CDO_
     - Climate Data Operators; command-line tools, including conservative re-mapping

   * - IDV_
     - more complete visualization

   * - NCO_
     - NetCDF Operators; command-line tools for pre- and post-processing

   * - NCL_
     - NCAR Command Language

   * - PyNGL_
     - Python version of NCL

The PISM authors use ``ncview`` and "``ncdump -h``" for quick visualization and metadata
examination. NCO has powerful command-line manipulation of NetCDF files, but requires some
learning. Another such command-line tool is CDO, but to use CDO on PISM files first run
the script ``nc2cdo.py``, from the ``util/`` PISM directory, on the file to fix the
metadata so that CDO will understand the mapping. Finally, Python scripts using the
``netcdf4-python`` package (see :ref:`sec-installation`) are often the best way to
non-trivially change a NetCDF file or make publishable figures from it. MATLAB also has
good NetCDF I/O capabilities.

See :numref:`tab-modelhierarchy` in section :ref:`sec-model-hierarchy` for an
overview on the data necessary for modeling. For more information on the format of input
files for PISM, see section :ref:`sec-initboot`.
