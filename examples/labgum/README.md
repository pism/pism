# laboratory validation example

Validation of isothermal SIA numerical model using laboratory
experiment.  The source of the set-up is

  R. Sayag and M. G. Worster, 2012. *Axisymmetric gravity currents of
  power-law fluids over a rigid horizontal surface*, J. Fluid Mech.
  doi:10.1017/jfm.2012.545, 1--9

At this point this example does not include the original laboratory data.

It is also unclear the exact dimension of the pipe supplying fluid.
The flow rate in the pipe is assumed to be constant across the pipe,
so that the `climatic_mass_balance` variable is constant in the pipe,
but some other profile could be built.

## Basic usage

Preprocess both translates `gumparams.cdl` --> `gumparams.nc` and it
creates `initgum.nc`:

    $ ./buildgum.py

Here is a run for 746 seconds on an 10 mm grid

    $ ./rungum.sh 4 51 &> out.lab51 &

When it is done, show the radial time series

    $ ./showradius.py -o r51.png ts_lab51.nc

Show the radial time series compared to data from Sayag (not currently public):

    $ ./showradius.py -o r51.png -d constantflux3.txt ts_lab51.nc

Results are much better on finer grids because the input pipe radius is
only 10 mm.

## Changing configuration constants

Edit `gumparams.cdl` and then rerun `buildgum.py` and then use `rungum.sh`.
