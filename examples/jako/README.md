Jakobshavn outlet glacier regional modeling example
=================

This directory contains all of the scripts needed to build a PISM regional model
of Jakobshavn Isbrae in the Greenland ice sheet.  The same strategy will work
for other outlet glaciers.  The base data is the SeaRISE 1km Greenland dataset
for the whole ice sheet, plus a 5km-resolution whole ice sheet spun-up state
called g5km_gridseq.nc.  That whole ice sheet state can be downloaded or
regenerated from examples/std-greenland/ scripts.

This example demonstrates the outlet-glacier-mode PISM regional mode `pism -regional`
and the drainage-basin-delineation tool `regional-tools`.

For more detail and complete instructions see the PISM User's Manual.
