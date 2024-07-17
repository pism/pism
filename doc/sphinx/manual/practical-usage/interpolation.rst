.. include:: ../../global.txt

.. default-role:: literal

.. _sec-interpolation:

Interpolation of input data
---------------------------

PISM includes a basic implementation of bilinear and trilinear interpolation. This code
assumes that input and output grids\ [#f1]_ use the same projection. (This is
interpolation in the projected *Cartesian* coordinate system.) It works best if input and
output grids use *similar* grid resolution. It suffers from signal aliasing when the input
grid resolution is significantly (e.g. 10 times) higher than the output grid resolution.\
[#f2]_

One way to avoid aliasing issues is by pre-processing PISM's inputs, e.g. *conservatively*
re-projecting and interpolating a data set from its original grid to the grid that will be
used by PISM. See CDO_ (specifically, `cdo remapcon`) and GDAL_ (`gdalwarp -r average` and
similar), for example. The downside of this approach is that one has to store and archive
pre-processed inputs.

It is possible to avoid these issues and perform re-projection and interpolation "on the
fly" by compiling PISM with PROJ_ and YAC_; see section :ref:`sec-install-yac`. (Note that
YAC performs interpolation *on a sphere* using longitudes and latitudes of cell centers
and corners.)

To use this code, make sure that

- the coordinate reference system used by PISM's internal grid was specified (see
  :ref:`sec-crs`)
- an input (often climate or ocean forcing) file contains projection information (again,
  see :ref:`sec-crs`).

PISM will automatically choose the interpolation method depending of the *interpolation
direction*:

- fine to coarse: first order conservative
- coarse to fine: distance-weighted sum of neighbors (similar to bilinear).

.. rubric:: Footnotes

.. [#f1] Here the "output" grid is the grid used to perform a simulation, i.e. PISM's
         internal grid.

.. [#f2] Aliasing artifacts in the bed topography field can significantly increase the
         "roughness" of the bed used by PISM, which unnecessarily increases the
         computational cost of a simulation.
