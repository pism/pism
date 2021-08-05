### A simplified setup showing that `-kill_icebergs` should always be "on".

We start with three blobs of ice on a flat bed.  As the run goes on,
the sea level rises and blobs start floating, starting with the
smallest one.

Once the whole blob is afloat the SSA solver will fail because of a
zero pivot.  The option `-kill_icebergs` removes such disconnected floating ice.

To run this example, do

    $ make

First look at variable `thk` in `flood.nc` to see initial geometry.
See evolving `mask` and `thk` variables in output files `ex_*.nc` to see
the effect of various options including `[none]`, `-cfbc`,
`-cfbc -kill_icebergs`.

Note that without the calving front boundary condition (`make no_cfbc`
which uses `[none]`), the strength extension "ties" floating blobs to the
biggest blob, which remains grounded.  In this (not-recommended) case the
SSA solver does not see a zero pivot.

As a result of this example and others, the option combination
`-cfbc -kill_icebergs` is generally recommended for marine ice sheet
simulations.
