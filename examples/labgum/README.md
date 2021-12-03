# laboratory validation example

This example is a validation of [PISM's](http://www.pism.io) isothermal
SIA numerical model using a laboratory experiment with a
[Xanthan gum](http://en.wikipedia.org/wiki/Xanthan_gum) suspension in water.
This fluid is more strongly shear-thinning than ice but it has nearly the same
density.  This example is documented in section 12.2 of the PISM User's Manual.

The source of the set-up is the "constant flux" experiment in

  R. Sayag and M. G. Worster, 2013. *Axisymmetric gravity currents of
  power-law fluids over a rigid horizontal surface*, J. Fluid Mech. 716,
  [doi:10.1017/jfm.2012.545](http://dx.doi.org/10.1017/jfm.2012.545).

See also

  R. Sayag, S. S. Pegler, and M. G. Worster, 2012. *Floating extensional flows*,
  Physics of Fluids 24 (9),
  [doi:10.1063/1.4747184](http://dx.doi.org/10.1063/1.4747184).

In the laboratory experiment, fluid is pushed through the bottom of a flat table
from a tube of radius 8 mm (R. Sayag, personal communication) at a mass rate of
about 3 g/s (Sayag & Worster, 2013).  The glaciological analog is an ice sheet
on a flat bed fed by positive surface mass balance in the vicinity of the dome,
but with zero surface mass balance everywhere else.

Sayag & Worster (2013) estimate n = 5.9 using regression of laboratory
measurements of the radius; see Figures 2(c) and 2(d) in the paper.  More
precisely, they compare radius data to a similarity solution of the thickness
evolution equation to infer the exponent n.

See `preprocess.py` for settings of various parameters which are appropriate to
this fluid problem.  The final total mass of fluid, about 1 kg, is about 18
orders of magnitude smaller than the mass of the Greenland ice sheet.)

The flow rate in the pipe is assumed to be constant across the pipe,
so the input "climate" has given `climatic_mass_balance` which is constant in
the pipe and zero outside the pipe.

## basic usage

The preprocessing stage builds a NetCDF file suitable for PISM bootstrapping.
It is built at exactly the run-time resolution in order to make the flux
into the center of the "ice" sheet have the correct value given the small size
of the positive mass flux area.  Also it creates `gumparams.nc` which contains
`pism_overrides` attributes, the special parameters for this experiment.

    $ ./preprocess.py

Now view `initlab52.nc`.  Only the `climatic_mass_balance` variable is
interesting.

Now we run for 746 model seconds (Sayag & Worster, 2013) on a 10 mm grid
(520 mm / 52 subintervals) using 4 processors:

    $ ./rungum.sh 4 52 &> out.lab52

This run generates `out.lab52` from `stderr` and `stdout`,
and it also generates diagnostic NetCDF files `ts_lab52.nc` and `ex_lab52.nc`.
It takes about 5 minutes on a 2013 laptop.

Results are better on finer grids because the input pipe radius is only 8 mm.
For example, this uses a 5 mm grid, and takes about an hour to run:

    $ ./preprocess.py -Mx 104 -o initlab104.nc
    $ ./rungum.sh 4 104 &> out.lab104

You can compare multiple runs to the experimental data on radius (R. Sayag,
personal communication):

    $ ./showradius.py -o foo.png -d constantflux3.txt ts_lab*.nc

## higher resolutions

For a 2.5 mm grid taking several hours to run, do:

    $ ./preprocess.py -Mx 208 -o initlab208.nc
    $ ./rungum.sh 4 208 &> out.lab208

More processors (NN) are recommended for a 1.0 mm grid:

    $ ./preprocess.py -Mx 520 -o initlab520.nc
    $ ./rungum.sh NN 520 &> out.lab520

To experiment with different configuration constants, edit and rerun
`preprocess.py`.

