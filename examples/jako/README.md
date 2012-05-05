Jakobshavn outlet glacier regional example
=================

This directory contains all of the scripts needed to build a PISM regional
of Jakobshavn Isbrae in the Greenland ice sheet.  The same strategy will
work for other outlet glaciers.  The base data is the SeaRISE 1km
Greenland dataset for the whole ice sheet plus an earlier lower-resolution
(5km) whole ice sheet spun-up state.  That whole ice sheet state can be
downloaded or regenerated from examples/searise-greenland/ scripts.

This example demonstrates the outlet-glacier-mode PISM executable `pismo`
and the drainage-basin-delineation tool `regional-tools`.


Get the drainage basin tool
----------

First get `regional-tools`, the drainage basin generator, using GIT.  Then
set up this python tool as directed in its `README.md`.  Then come back here and
link the tool we need.  Here is the quick summary of these actions:

    $ cd ~                                 # whatever location you want for it
    $ git clone https://github.com/pism/regional-tools.git
    $ cd regional-tools/
    $ python setup.py install              # or add "--user" to build locally
    $ cd ~/pism/examples/jako/             # come back to current directory
    $ ln -s ~/regional-tools/pism_regional.py .      # symbolic link to tool

To see a description of the drainage basin tool itself, and a bit on how it
works, see https://github.com/pism/regional-tools.


Preprocess fine grid geometry data and coarse climate data
----------

Next we use a script that downloads and cleans the SeaRISE data, an 80 Mb file
called `Greenland1km.nc`.

    $ ./preprocess.sh

If the file is already present then no download occurs, but preprocessing
proceeds in any case.

When it is done, a file `gr1km.nc` is created.  We can look at what fields
it contains:

    $ ncdump -h gr1km.nc                   # view metadata
    $ ncview gr1km.nc                      # view fields

Also we download the SeaRISE 5km data set `Greenland_5km_v1.1.nc` because it
contains the surface mass balance model result from RACMO.  This field is not present
in the SeaRISE 1km data set.  If you have already run the example in the first
section of the PISM User's Manual, i.e. examples/searise-greenland/ then you
already have this file and you can link to it to avoid downloading:

    $ ln -s ../searise-greenland/Greenland_5km_v1.1.nc

But another preprocess stage is still needed to fix the metadata, so do:

    $ ./getsmb.sh


Get coarse grid, whole ice sheet spinup result
-----------

Finally, we also get coarse grid model results for the whole ice sheet.  These are 
not data but come from a prior PISM run.  They provide the boundary conditions
and thermodynamical initial condition for 
the regional flow model we are building.  The next script downloads a
precomputed whole-ice-sheet result from a PISM SeaRISE-Greenland run, a 1.3Gb
file called `g5km_0_ftt.nc`.  (No download occurs if a file by that name is
already present.)

    $ ./getcoarse.sh

A file `g5km_bc.nc` is created with the fields we actually need.  In particular
there are thermodynamical spun-up variables (`enthalpy`,`bmelt`,`bwat`) and
boundary values for the sliding velocity in the regional model (`u_ssa_bc`,
`v_ssa_bc`).


Creating a region for us to model
-------------

We are going to extract a drainage basin mask from the fine grid surface
elevation data.  The outline of the flow into the outlet glacier we want to model
can be determined by the gradient flow from the surface elevation.
Within this masked drainage basin we will apply all physics in the PISM model---
see below---but outside this basin we will apply simplified models
or use precomputed whole ice sheet results as boundary conditions.

So we use the script `pism_regional.py` from `regional-tools` (see above)
on `gr1km.nc` (created above):

    $ python pism_regional.py

Running the script opens a graphical user interface
(GUI) which allows you to select a small rectangle around the flux gate or
fjord which is the terminus of the outlet glacier you want to model.

**To use the GUI:**  Select `gr1km.nc` to open.  Once the topographic map appears
in the Figure 1 window, select the button `Select terminus rectangle`.
Now use the mouse to select a small rectangle around the Jakobshavn
terminus (calving front).  (<em>At this stage you can choose any other
terminus/calving front and put a small rectangle around it!</em>)
Once you have a highlighted rectangle, select a `border width` of at
least 50 cells.  (<em>This suggestion is somewhat Jakobshavn-specific.
The intention is to have an ice-free western boundary on or region.</em>)

Then click `Compute the drainage basin mask`.  Because this is a large data set
there will be some delay.  (<em>A parallel computation is done.</em>)

Finally click `Save the drainage basin mask` and save with your
preferred name; we will assume here that it is `jakomask.nc`.  Then quit.  You
can look at the result with `ncview`:

    $ ncview -minmax all jakomask.nc


Re-creating the region mask if needed
-------------

For repeatability the rest of this **README.md** assumes a particular
choice of region.  That is, you can skip the GUI usage above and run

    $ python pism_regional.py -i gr1km.nc -o jakomask.nc -x 360,382 -y 1135,1176 -b 50

This also generates the file `jakomask.nc` used in the rest of the script.

Such a step is exactly what is needed to have more precise control over the
identification of a terminus rectangle.  You may need to re-create the region
with slight modifications, for instance.  For more options:

    $ python pism_regional.py --help


Cut out our region from each whole ice sheet data set
----------

We still need to "cut out" the modeled region we want from the whole sheet data
sets `gr1km.nc` and `g5km_bc.nc`.  (The coarse grid climate data `g5km_climate.nc`
does not need this action because PISM's coupling code can already handle all
needed interpolation/subsampling for such climate data.)

You may have noticed that the feedback from the previous stage includes a
cutout command which appears as a global attribute of jakomask.nc.  Get it
and use it:

    $ ncdump -h jakomask.nc |grep cutout

Cut and paste the command, which applies to both the 1km Greenland data file
and the mask file.  For example it will look like:

    $ ncks -d x,299,918 -d y,970,1394 gr1km.nc jako.nc
    $ ncks -A -d x,299,918 -d y,970,1394 jakomask.nc jako.nc   # note -A for append!

Now look at `jako.nc`:

    $ ncview jako.nc

This file is the full geometry data ready for a regional model.  The field
`ftt_mask` has an identified drainage basin, outside of which we
will not use a full model, but use simplified time-independent boundary
conditions instead.  Specifically, outside of the `ftt_mask` area we will
essentially keep the initial thickness, while inside the `ftt_mask` area all
fields will evolve normally.

To run PISM we will need to know the size of the 1km grid in jako.nc.  Do this:

    $ ncdump -h jako.nc |head -n 4
    netcdf jako {
    dimensions:
	  y = 425 ;
	  x = 620 ;

That is, our region covers a 620 km by 425 km region of Greenland.  A 1km
resolution, century scale model run is achievable on a desktop machine, and 
that is our goal below.  (<em>PISM can do 1km runs for the whole ice sheet,
but you do need a supercomputer for that!</em>)


Running the model on a 3km grid
-----------

The whole ice sheet model fields in `g5km_bc.nc` may or may not, depending on
modeler intent, be spun-up adequately for the purposes of the regional model.
That is, the intent may be to study equilibrium states with model settings
special to the region.  Here we assume that is the case.  Thus we do more
spin-up and we attempt to get an equilibrium 1km model.

In the first stage we therefore work on an intermediate grid between the 1km
goal and the 5km whole ice sheet state.  Quick calculations (round up on
`620/3 + 1` and `425/3 + 1`) suggest `-Mx 208 -Mx 143` as a 3km grid for PISM.

So now we do a basic run using 4 MPI processes:

    $ ./spinup.sh 4 208 143 >> out.spin3km &

Please read the script, and/or watch the run while it gets going:

    $ less spinup.sh
    $ less out.spin3km

**Comments:** There are two runs, first a 1yr run with just SIA velocities,
which sligin  FIXME:  maybe that's stupid

FIXME FROM HERE

The runs get Dirichlet BCs from `g5km_bc.nc`.  Specifically there
is option `-no_model_strip 10` asking `pismo` to put a 10km strip around the
computational domain wherein the thermodynamical spun-up variables
(`enthalpy`,`bmelt`,`bwat`) are held fixed and used as boundary conditions for
the conservation of energy model.  The boundary values for the sliding velocity
(`u_ssa_bc`,`v_ssa_bc`) are also held fixed in the same strip.

(<em>An alternative is to have the enthalpy not spun-up at all in the strip,
which would happen if the first run 
which then advects inward as an unrealistic temp.  The point for the SSA bc
is that the alternative is to have zero velocity there while the velocity 
tangent to the edge of the no_model strip, i.e. in westerly direction, is
significantly nonzero(?).



The model on the 1km grid
-----------

FIXME Thus here is a native 1km one year run: FIXME

