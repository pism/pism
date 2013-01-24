EISMINT II helper script
=========

The PISM User's Manual documents the EISMINT II experiments.

The only thing here is a helper script `runexp.sh`.

Run as

    $ ./runexp.sh NN EXP MM [INPUTFILE]

for `NN` procs, experiment `EXP` and an `MM x MM x MM` grid in three dimensions.
The fourth (optional) argument specifies an input file.  For example,

    $ ./runexp.sh 4 A 61

is the standard EISMINT II experiment A run on the standard grid.  A final-state
output file `eisIIA61.nc` and diagnostic files `ts_eisIIA61.nc` and
`ex_eisIIA61.nc` will be produced.

Do this for a dry run

    $ PISM_DO=echo ./runexp.sh 4 A 61

Specify a NetCDF file name as the fourth argument for a run starting from an
input file.  For example the standard EISMINT II experiment B run should start
with the final state of experiment A:

    $ ./runexp.sh 4 B 61 eisIIA61.nc
