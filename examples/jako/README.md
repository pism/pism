Jakobshavn outlet glacier regional example
=================

This directory contains all of the scripts needed to build a PISM regional
of Jakobshavn Isbrae in Greenland.  The same strategy should work for most
other outlet glaciers in that ice sheet.  The base data is the SeaRISE 1km
Greenland dataset for the whole ice sheet.

This example demonstrates the outlet-glacier-mode PISM executable `pismo`
and the drainage-basin-delineation tool 

Preprocessing
----------

In order to run `pismo` some preprocessing must be performed.  A number
of python tools are needed.

First get `regional-tools`, the drainage basin generator, by `git clone`.  Then
set up this python tool as directed in its `README.md`.  Then come back here and
link the tool we need.  Here is the quick summary:

    $ cd ~                                 # location where you want it
    $ git clone https://bueler@github.com/pism/regional-tools.git
    $ cd regional-tools/
    $ python setup.py install              # or add "--user" to build locally
    $ cd ~/pism/examples/jako/             # directory with current README
    $ ln -s ~/regional-tools/pism_regional.py .

Next we use a script that downloads and cleans the SeaRISE data, an 80 Mb file called
`Greenland1km.nc`.

    $ ./preprocess.sh

A file `gr1km.nc` is created.  We can look at what fields it contains:

    $ ncdump -h gr1km.nc                   # view metadata
    $ ncview gr1km.nc                      # view fields

We also download coarse grid information for the whole ice sheet.  These are 
model results already, which provide the boundary conditions and climate for 
the regional flow model we are building.  The next script downloads a
precomputed whole-ice-sheet result from a PISM SeaRISE-Greenland run, a 1.3Gb
file called `g5km_0_ftt.nc`.  A file `g5km_bc.nc` is created with the fields
we actually need.

Also we download the SeaRISE 5km data set `Greenland_5km_v1.1.nc` because it
contains the surface mass balance model result from RACMO, which is not present
in the SeaRISE 1km data set.  If you have already run the example in the first
section of the PISM User's Manual then you already have this file and you can
link to it to avoid downloading.  Thus

    $ ln -s ../searise-greenland/Greenland_5km_v1.1.nc
    $ ./getcoarse.sh


Creating a regional model
-------------

First we extract a drainage basin mask from the surface elevation data in
`gr1km.nc`.  The idea is that the outline of the flow we want to model can be
determined by the gradient flow from the surface elevation.  Then
within the masked drainage basin we will apply all physics in the PISM model
which we run below.  But outside the masked basin we apply simplified models
or precomputed whole ice sheet results as boundary conditions.

So we use `pism_regional.py` on `gr1km.nc`.  This opens a graphical user interface
(GUI) which allows you to select a small rectangle around the flux gate or
fjord which is the terminus of the outlet glacier you want to model.

    $ python pism_regional.py

To use the GUI:  Select `gr1km.nc` to open.  Once the topographic map appears
in the Figure 1 window, select the button `Select terminus rectangle`.
Now use the mouse to select a small rectangle around the Jakobshavn
terminus, for example.  Once it is a highlighted rectangle, click
`Compute the drainage basin mask`.  Because this is a large data set
there will be some delay.  (Note that a parallel computation is done.)

Finally click `Save the drainage basin mask` and save with your
preferred name; we will assume here that it is `jakomask.nc`.

  FIXME:
    *  selecting rectangle slow/fragile?
    *  feedback on completion of steps should appear in GUI not in terminal?
    *  cutout command which appears at terminal can get lost?
    *  I have to hit "quit" twice to make it so

We still need to "cut out" the region we want, including the whole
drainage basin.
The feedback from the above process includes a cutout command which
appears as a global attribute of jakomask.nc.  Get it and use it:

    $ ncdump -h jakomask.nc |grep cutout

Cut and paste the command, which applies to both the 1km Greenland data file
and the mask file.  For example it will look like:

    $ ncks -d x,340,873 -d y,1015,1349 gr1km.nc jako.nc
    $ ncks -A -d x,340,873 -d y,1015,1349 jakomask.nc jako.nc   # append

Now look at `jako.nc`.  This file is ready for a regional model.  The field
`ftt_mask` has an identified drainage basin, outside of which we
will not use a full model, but use simplified time-independent boundary
conditions instead.  Specifically, outside of the `ftt_mask` area we will
essentially keep the initial thickness, while inside the `ftt_mask` area all
fields will evolve normally.

FIXME What is the grid in jako.nc?

    $ ncdump -h jako.nc |head -n 5
    netcdf jako {
    dimensions:
	  y = 335 ;
	  x = 534 ;
    variables:

FIXME Thus here is a native 1km one year run:  [extract from spinup3km.sh]

Running the model
-----------

To make real progress on a spinup we divide the grid dimensions by 3 to get a
3km grid.  On this grid we do a basic run using 4 MPI processes:

    $ ./spinup.sh 4 >> out.spin3km &

Comments: The run gets Dirichlet BCs from a saved whole ice sheet run (g5km_0_ftt.nc
saved at pismdocs.org) for bmelt, enthalpy, SSA.  The point for the second of these is that
the only alternative is to have the enthalpy not spun-up at all in the no_model_strip,
which then advects inward as an unrealistic temp.  The point for the SSA bc
is that the alternative is to have zero velocity there while the velocity 
tangent to the edge of the no_model strip, i.e. in westerly direction, is
significantly nonzero(?).

