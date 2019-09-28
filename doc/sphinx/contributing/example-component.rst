.. include:: ../global.txt

.. default-role:: literal

.. _sec-writing-pism-components:

Adding a new PISM component
===========================

This section describes the implementation of a simple ocean model component.

Most of PISM's sub-models are derived from the abstract class `Component`. This class
contains data members used in most sub-models and defines a common interface. Each
sub-model has to be able to do the following.

#. Initialize using an input file.

#. Save its model state to an output file.

   To avoid a performance penalty when saving to a NetCDF-3 file this is performed in two
   steps:

   #. Define NetCDF variables, creating the file structure.
   #. Write data to the file.

#. Time-stepping sub-models need to provide the maximum allowable time step length at a
   given model time.

#. Provide lists of spatially-variable and scalar (time-series) diagnostics.


FIXME: add more content.

FIXME: write about the importance of getting initialization right.
