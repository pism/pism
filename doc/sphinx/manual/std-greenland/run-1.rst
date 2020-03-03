.. include:: ../../global.txt

.. _sec-runscript:

First run
---------

Like many Unix programs, PISM allows a lot of command-line options. In fact, because the
variety of allowed ice sheet, shelf, and glacier configurations, and included sub-models,
is so large, command-line options are covered in sections :ref:`sec-initboot` through
:ref:`sec-practical-usage` of this manual.\ [#]_ In practice one often builds scripts to run
PISM with the correct options, which is what we show here. The script we use is
"``spinup.sh``" in the ``examples/std-greenland/`` subdirectory of ``pism/``.

Note that initializing ice sheets, usually called "spin-up", can be done by computing
approximate steady states with constant boundary data, or, in some cases, by integrating
paleo-climatic and long-time-scale information, also applied at the ice sheet boundary, to
build a model for the present state of the ice sheet. Both of these possibilities are
illustrated in the ``spinup.sh`` script. The spin-up stage of using an ice sheet model may
actually require more processor-hours than follow-on "experiment" or "forecast" stages.

To see what can be done with the script, read the usage message it produces:

.. code-block:: none

   ./spinup.sh

The simplest spin-up approach is to use a "constant-climate" model. We take this approach
first. To see a more detailed view of the PISM command for the first run, do:

.. literalinclude:: scripts/run-1-echo.sh
   :language: bash
   :lines: 3-

Setting the environment variable ``PISM_DO`` in this way tells ``spinup.sh`` just to print
out the commands it is about to run instead of executing them. The "proposed" run looks
like this:

.. literalinclude:: scripts/run-1-command.sh
   :language: bash
   :name: firstcommand
   :lines: 3-

Let's briefly deconstruct this run.

At the front is "``mpiexec -n 4 pismr``". This means that the PISM executable ``pismr`` is
run in parallel using four processes (usually one per CPU core) under the `Message Passing
Interface <MPI_>`_. Though we are assuming you have a workstation or laptop with at least
4 cores, this example will work with 1 to about 50 processors, with reasonably good
scaling in speed. Scaling can be good with more processors if we run at higher spatial
resolution :cite:`BBssasliding`, :cite:`DickensMorey2013`. The executable name "``pismr``"
stands for the standard "run" mode of PISM (in contrast to specialized modes described
later in sections :ref:`sec-verif` and :ref:`sec-simp`).

Next, the proposed run uses option ``-bootstrap`` to start the run by "bootstrapping."
This term describes the creation, by heuristics and highly-simplified models, of the
mathematical initial conditions required for a deterministic, time-dependent ice dynamics
model. Then the options describe a `76 \times 141` point grid in the horizontal,
which gives 20 km grid spacing in both directions. Then there are choices about the
vertical extent and resolution of the computational grid; more on those later. After that
we see a description of the time axis, with a start and end time given: "``-ys -10000 -ye
0``".

Then we get the instructions that tell PISM to read the upper surface boundary conditions
(i.e. climate) from a file: "``-surface given -surface_given_file
pism_Greenland_5km_v1.1.nc``". For more on these choices, see subsection
:ref:`sec-climate-inputs`, and also the :ref:`Climate Forcing Manual
<sec-climate-forcing>`.

Then there are a couple of options related to ice dynamics. First is a minimal "calving"
model which removes ice leaving the area covered by ice (or ice-free with the bed above
the sea level) according to the input file ("``-front_retreat_file``"); see section
:ref:`sec-calving` for this and other calving options). Then there is a setting for
enhanced ice softness ("``-sia_e 3.0``"). See section :ref:`sec-rheology` for more on this
enhancement parameter, which we also return to later in section :ref:`sec-paramstudy`.

Then there are longish options describing the fields we want as output, including scalar
time series ("``-ts_file ts_g20km_10ka.nc -ts_times -10000:yearly:0``"; see section
:ref:`sec-practical-usage`) and space-dependent fields ("``-extra_file ...``"; again see
section :ref:`sec-practical-usage`), and finally the named output file ("``-o
g20km_10ka.nc``").

Note that the modeling choices here are reasonable, but they are not the only way to do
it! The user is encouraged to experiment; that is the point of a model.

Now let's actually get the run going:

.. literalinclude:: scripts/run-1.sh
   :lines: 3-
   :language: bash

The terminating "``&``", which is optional, asks the system to run the command in the
background so we can keep working in the current shell. Because we have re-directed the
text output ("``&> out.g20km_10ka``"), PISM will show what it is doing in the text file
``out.g20km_10ka``. Using ``less`` is a good way to watch such a growing text-output file.
This run should take around 20 minutes.

.. rubric:: Footnotes

.. [#] Moreover, every configuration parameter can be set using a command-line option with
       the same name.
