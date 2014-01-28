# laboratory validation example

This example is a validation of [PISM's](http://www.pism-docs.org) isothermal
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
  Physics of Fluids 24 (9), [doi:10.1063/1.4747184](http://dx.doi.org/10.1063/1.4747184).

They push the fluid through the bottom of a flat table in a tube of radius
8 mm (R. Sayag, personal communication) at a mass rate of about 3 g / s
(Sayag & Worster, 2013).  The glaciological analog is an ice sheet on a flat bed
fed by positive surface mass balance in the vicinity of the dome, but with zero
surface mass balance everywhere else.

They estimate n = 5.9 using regression of laboratory measurements of the radius;
see Figures 2(c) and 2(d) in Sayag & Worster (2013).  More precisely, they use
radius data compared to a similarity solution of the geometry evolution
equation to infer the exponent n.  (Compare n = 3 for ice.)

See `gumparams.cdl` for settings of various parameters which suit this fluid
problem of total mass ~ 1 kg.  (Compare ~ 10^18 kg for the Greenland ice sheet.)

The flow rate in the pipe is assumed to be constant across the pipe,
so the input "climate" has given `climatic_mass_balance` which is constant in
the pipe and zero outside the pipe.

## basic usage

The preprocessing stage builds a NetCDF file suitable for PISM bootstrapping.
It is built at exactly the run-time resolution in order to make the flux
into the center of the "ice" sheet have the correct value given the small size
of the positive mass flux area.  Also it converts `gumparams.cdl` to
`gumparams.nc`

    $ ./preprocess.py

Now view `initlab53.nc`; only the `climatic_mass_balance` variable is
interesting.  Now we run for 746 model seconds (compare Sayag & Worster (2013))
on an 10 mm grid (520 mm / 52 subintervals) using 4 processors:

    $ ./rungum.sh 4 53 &> out.lab53 &

This run generates `out.lab53` from `stdout`
and it also generates diagnostic NetCDF files `ts_lab53.nc` and `ex_lab53.nc`.
It took about 5 minutes on a 2013 laptop.  When it is done, you can show the
modeled radius compared to the experimental data (R. Sayag, personal
communication):

    $ ./showradius.py -o r53.png -d constantflux3.txt ts_lab53.nc

Results are better on finer grids because the input pipe radius is only 8 mm.
For example, this uses a 5 mm grid, and takes about an hour to run:

    $ ./preprocess.py -Mx 105 -o initlab105.nc
    $ ./rungum.sh 4 105 &> out.lab105 &

while the next uses a 2.5 mm grid, and takes several hours to run:

    $ ./preprocess.py -Mx 209 -o initlab209.nc
    $ ./rungum.sh 4 209 &> out.lab209 &

Note you can compare multiple runs to the data in one figure:

    $ ./showradius.py -o foo.png -d constantflux3.txt ts_lab*.nc

To experiment with different configuration constants, edit `gumparams.cdl` and
then rerun `preprocess.sh`.
