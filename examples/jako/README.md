Jakobshavn outlet glacier regional modeling example
=================

This directory contains all of the scripts needed to build a PISM regional model
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
    $ python setup.py install              # add "--user" to build locally
    $ cd ~/pism/examples/jako/             # come back to current directory
    $ ln -s ~/regional-tools/pism_regional.py .      # symbolic link to tool

To see a description of the drainage basin tool itself, and a bit on how it
works, see https://github.com/pism/regional-tools.


Preprocess higher-resolution geometry data
----------

Next we use a script that downloads and cleans the fine-grid SeaRISE data,
an 80 Mb file called `Greenland1km.nc`.

    $ ./preprocess.sh

If the file is already present then no download occurs, but preprocessing
proceeds in any case.

When it is done, a file `gr1km.nc` is created.  We can look at what fields
it contains:

    $ ncdump -h gr1km.nc                   # view metadata
    $ ncview gr1km.nc                      # view fields


Preprocess lower-resolution climate data
----------

Also we download the SeaRISE 5km data set `Greenland_5km_v1.1.nc` because it
contains the surface mass balance model result from RACMO.  This field is not present
in the SeaRISE 1km data set above.  If you have already run the example in the first
section of the PISM User's Manual, i.e. in `examples/searise-greenland/`, then you
already have this file and you can link to it to avoid downloading:

    $ ln -s ../searise-greenland/Greenland_5km_v1.1.nc

In any case, run this additional preprocess stage to fix the metadata:

    $ ./getsmb.sh

A file `g5km_climate.nc` will appear, and can be examined in the usual ways.


Get whole ice sheet spinup result
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


Extract region for modeling
-------------

We are going to extract a "drainage basin mask" from the fine grid surface
elevation data.  The outline of the flow into the outlet glacier we want to model
can be determined by the gradient flow from the surface elevation.
`pism_regional.py` will identify the upstream area, the drainage basin in the
sense of the surface gradient flow, into a chosen "terminus rectangle".
Within this masked drainage basin we will apply all physics in the PISM model
but outside this basin we will apply simplified models
or use precomputed whole ice sheet results as boundary conditions.

So we use the script `pism_regional.py` from `regional-tools` (see above)
on `gr1km.nc` (created above):

    $ python pism_regional.py

Running the script opens a graphical user interface
(GUI) which allows you to select a small rectangle around the flux gate or
fjord which is the terminus of the outlet glacier you want to model.

**To use the GUI:**  Select `gr1km.nc` to open.  Once the topographic map appears
in the Figure 1 window, you may zoom enough to see the general outlet glacier
area.  Then select the button `Select terminus rectangle`.
Now use the mouse to select a small rectangle around the Jakobshavn
terminus (calving front).  (<em>At this stage you can choose any other
terminus/calving front and put a small rectangle around it!</em>)
Once you have a highlighted rectangle, select a `border width` of at
least 50 cells.  (<em>This suggestion is somewhat Jakobshavn-specific.
The intention is to have an ice-free western boundary on the computational
domain for our modeled region.</em>)

Then click `Compute the drainage basin mask`.  Because this is a large data set
there will be some delay.  (<em>A parallel computation is done.</em>)

Finally click `Save the drainage basin mask` and save with your
preferred name; we will assume here that it is `jakomask.nc`.  Then quit.  You
can look at the result with `ncview`:

    $ ncview -minmax all jakomask.nc

**Re-creating the region mask if needed:** For repeatability the rest
of this `README.md` assumes a particular
choice of region.  That is, you can skip the GUI usage above and run

    $ python pism_regional.py -i gr1km.nc -o jakomask.nc -x 360,382 -y 1135,1176 -b 50

The options `-x A,B -y C,D` identify the grid indices of the corners of the
terminus rectangle.  Thus this command also generates the file `jakomask.nc`
used in the rest of the script.

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

You may have noticed that the text output from `pism_regional.py` includes a
cutout command which appears as a global attribute of jakomask.nc.  Get it this way:

    $ ncdump -h jakomask.nc |grep cutout

Copy and run the command that appears in the string `cutout`, something like

    $ ncks -d x,299,918 -d y,970,1394 gr1km.nc jako.nc

This command is applied to both the 1km Greenland data file and the
mask file, so modify the command to work on the mask also, and note option
`-A` for append:

    $ ncks -A -d x,299,918 -d y,970,1394 jakomask.nc jako.nc

Now look at `jako.nc`:

    $ ncview -minmax all jako.nc

This file is the full geometry data ready for a regional model.  The field
`ftt_mask` has an identified drainage basin, outside of which we
will not use a full model, but use simplified time-independent boundary
conditions instead.  Specifically, outside of the `ftt_mask` area we will
essentially keep the initial thickness, while inside the `ftt_mask` area all
fields will evolve normally.

To run PISM we will need to know the size of the 1km grid in jako.nc.  Do this:

    $ ncdump -h jako.nc |head
    netcdf jako {
    dimensions:
	  y = 425 ;
	  x = 620 ;
	...

The grid has spacing of 1 km, so our region is a 620 km by 425 km rectangle.


Spinning-up the regional model on a 5km grid
-----------

A 2km resolution, century-scale model run on this Jakobshavn region is
achievable on a desktop machine, with a bit of patience, and that is our goal below.
A 5km resolution run, which is the resolution of the 5km whole ice sheet state
which was computed on a supercomputer, is definitely achievable on a
desktop or laptop, and that is what we do first.

(<em>PISM can do 1km runs for the whole Greenland ice sheet, as shown in this
[`pism-docs.org` news item](http://www.pism-docs.org/wiki/doku.php?id=news:first1km)
but you certainly need a supercomputer for that!</em>)

The whole ice sheet model fields in `g5km_bc.nc`, which is restricted to the 
regional grid as the initial state of the regional model, may or may not,
depending on modeler intent, be spun-up adequately for the purposes of the
regional model.  (For instance, the intention may be to study equilibrium states
with model settings special to the region.)  Here we assume that more
spin-up is needed.  We will get first an equilibrium 5km model, and then do a
century run of a 2km model based on the equilibriated 5km result.
(Determining "equilibrium" requires a decision, of course.  The standard here
is that the ice volume in the region should not change more than one percent
in 100 model years.  See `ivol` in `ts_spunjako_0.nc` below.)

Quick calculations for the 5km grid, e.g. by calculations
`620/5 + 1` and `425/5 + 1`, suggest `-Mx 125 -My 86`.
So now we do a basic run using 2 MPI processes:

    $ ./spinup.sh 2 125 86 >> out.spin5km &

Please read the script, and watch the run, while it runs:

    $ less spinup.sh
    $ less out.spin5km

The run takes about FIXME (approx 4?) processor-hours.

**Comments on the run:**
The run uses `-boot_file` on `jako.nc`.  A modestly-fine vertical
grid (20 m spacing) is chosen.

There is option `-no_model_strip 10` asking
`pismo` to put a 10km strip around edge of the computational domain wherein
some variables will be held fixed.  Specifically the thermodynamical spun-up
variables `enthalpy`,`bmelt`,`bwat` from `g5km_bc.nc` are held fixed and used
as boundary conditions for the conservation of energy model in the strip.
Additionally, Dirichlet boundary conditions `u_ssa_bc`,`v_ssa_bc` are read
for the sliding stress balance (the SSA) from the same file.

(<em>An alternative is to have the enthalpy and other thermodynamical variables
not spun-up at all in the strip, which would happen if options `-regrid_...`
were not used.  However, the resulting not-very-realistic ice temperatures and
softness/hardness is advected inward. An alternative for the SSA boundary
conditions is to have zero velocity in the `no_model_strip`, but the velocity
tangent to the north and south edges of the strip is significantly nonzero
in fact.  Thus, generally the regridding techniques used here are recommended
for regional modeling.</em>)

The calving front of the glacier is handled by the following options combination:
`-ocean_kill FILENAME -cfbc -kill_icebergs`.  This choice uses the present-day
ice extent to determine the location of the calving front, but it then applies
the PIK mechanisms for the stress boundary condition at the calving front
(`-cfbc`) and it uses `-kill_icebergs` to eliminate any stray floating pieces
of ice for which stress balances are indeterminant (ill-posed).


Century run on a 2km grid
-----------

FIXME

    $ ./century.sh 8 311 212 spunjako_0.nc >> out.2km_100a &

