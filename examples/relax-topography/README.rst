Example: relax Greenland's basal topograpy
======================

This minimal example shows how we can calculate the isostatically-relaxed
basal topography for an ice-free Greenland.  The same method can be
applied to other areas.  It uses the flux correction method
"force-to-thickness" with a target thickness of zero to remove all ice
combined with the Lingle-Clark bed deformation model to calculate
isostatic uplift.

How to run it
-------------------------

First, download the SeaRISE data set running

``$ ./preprocess.sh``

or, alternatively, provide a data set of your choice.  For example, you can
make a link to the SeaRISE data set downloaded into
``examples/searise-greenland/``, if the example in Chapter 1 of the
User's Manual has already run.

Then run

``$ ./run-relax.sh N M``

where N is the number of processors and M defines the grid resolution
(0: 20km [default], 1: 10km, 2: 5km, 3: 2.5km, 4: 2km, 5: 1km). To run the
example on 16 cores and 5km horizontal grid resolution, in the background,
and store ``stdout`` output in a text file, type

``$ ./run-relax.sh 16 2 &> out.relax &``

First, the model is run for 1,000 years using flux-correction to
remove (most) of the ice. The second simulation continues from the
first, but now with the option -no_mass, for another 50,000 years
(Note: this number was defined ad-hoc, and may or may not be appropriate).
