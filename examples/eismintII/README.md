EISMINT II helper script
=========

The PISM User's Manual documents the EISMINT II experiments.  The only thing
here is a helper script `runexp.sh`.

Run as

    $ ./runexp.sh NN EXP MHOR MVER DUR

for `NN` procs, experiment `EXP`, an `MHOR x MHOR x MVER` grid in three
dimensions, and a run of `DUR` years.  All options are optional; see the script
for default settings.  Alternatively, run as

    $ ./runexp.sh NN EXP X X DUR [INPUTFILE]

The fifth argument specifies an input file; in this case the third and fourth
options `X X` are ignored.

For example,

    $ ./runexp.sh 4 A 61 121 1e4

is a shortened EISMINT II experiment A run on the standard horizontal grid but a
refined vertical grid.  A final-state output file `eisIIA61.nc` and diagnostic
files `ts_eisIIA61.nc` and `ex_eisIIA61.nc` will be produced.

Do this for a dry run

    $ PISM_DO=echo ./runexp.sh 4 A 61 121 1e4

Specify a NetCDF file name as the fourth argument for a run starting from an
input file.  For example this shortened EISMINT II experiment B run starts
with the final state of experiment A:

    $ ./runexp.sh 4 B X X 1e4 eisIIA61.nc
