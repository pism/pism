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

Here is a run for 300 seconds on an 8 mm grid

    $ ./rungum.sh 4 51 &> out.lab51 &


## Changing configuration constants

Edit `gumparams.cdl` and then rerun `buildgum.py` and then use `rungum.sh`.
