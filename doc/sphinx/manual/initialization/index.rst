.. include:: ../../global.txt

.. _sec-initboot:

Initialization and bootstrapping
================================

There are three ways to start PISM:

- option :opt:`-i` reads a previously-saved "complete" PISM model state from a NetCDF file,
  or
- option :opt:`-i` :opt:`-bootstrap` reads an "incomplete" NetCDF file and uses heuristics
  to fill in needed fields, or
- ``pismr -eisII ...`` and the ``pismv`` executable are used to initialize
  simplified-geometry experiments and verification tests from formulas in the source code,
  and thus no input file is required.

One of the first two choices is required when using the executable ``pismr``. Modeling
usually starts with the ``-i input.nc -bootstrap`` because real ice sheet observations are
never complete initial conditions. Runs with multiple stages often use the :opt:`-i`
option after the first stage.

.. toctree::

   restarting.rst

   bootstrapping.rst
