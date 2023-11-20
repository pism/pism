Ross flow model example
=================

This example demonstrates PISM modeling of ice shelves.  This model can be thought of as a validation case (see subsection 12.2 of the PISM User's Manual) or as a regional model.

We use a 5km data set for the entire Antarctic ice sheet.  Thus, based on the steps here, one can build models of any Antarctic ice shelf, either diagnostic (ice thickness prescribed, velocity boundary conditions prescribed, but computing flow velocity of the floating ice) or prognostic (ice thickness evolution).

These publications cover PISM applications to Antarctic ice shelves:

* T. Albrecht, A. Levermann (2012).  *Fracture field for large-scale ice dynamics*.  Journal of Glaciology 58 (207), 165--176.
* A. Levermann, T. Albrecht, R. Winkelmann, M. A. Martin, M. Haseloff, I. Joughin (2012) [*Kinematic first-order calving law implies potential for abrupt ice-shelf retreat*](http://www.the-cryosphere.net/6/273/2012/tc-6-273-2012.html).  The Cryosphere 6, 273--286.
* M. A. Martin, R. Winkelmann, M. Haseloff, T. Albrecht, E. Bueler, C. Khroulev, A. Levermann (2011).  [*The Potsdam Parallel Ice Sheet Model (PISM-PIK)--Part 2: Dynamic equilibrium simulation of the Antarctic ice sheet*](http://www.the-cryosphere.net/5/727/2011/tc-5-727-2011.html). The Cryosphere 5, 727--740.

Input datasets
==========

We use two public datasets which are 16 Mb and 96 Mb, respectively:

* [*An improved Antarctic dataset for high resolution numerical ice sheet models (ALBMAP v1)*](http://doi.pangaea.de/10.1594/PANGAEA.734145)  by A. M. Le Brocq, A. J. Payne, and A. Vieli
* [*MEaSUREs InSAR-Based Antarctica Velocity Map*](http://nsidc.org/data/nsidc-0484.html) by Rignot, E., J. Mouginot, and B. Scheuchl. 2011,

From these data sets, the following fields from ALBMAP are used:

* ice thickness (equivalently: the upper/lower surface elevation pair),
* ice surface temperature,
* ice surface mass balance,
* bedrock topography, and

From the MEaSUREs data set these fields are used:

* u and v components of the surface ice velocity for grounded ice

These velocity components become the horizontal boundary conditions for the floating ice shelf.  The values from the interior of the ice shelf are used only by `plot.py` (below) to evaluate the model velocities.


Running the examples
=================

First run

    $ ./preprocess.py

This will download the data sets `ALBMAPv1.nc.zip` and `antarctica_ice_velocity.nc.gz` if they are not already downloaded.  Then it uses NCO and CDO tools to fix metadata.

Now, depending on your intent you will choose among the three subdirectories:

* `diagnostic/`: Use PISM's SSA stress balance to compute flow velocity from geometry, ice hardness (a function of temperature), and observed velocities at the grounding line as boundary conditions.  This case is documented in subsection 12.2 of the PISM User's Manual.  See `README.m` and `run_diag.sh` in the subdirectory.

* `prognostic/`: Use PISM's SSA stress balance to do a time-stepping run, which computes evolcing flow velocity and evolving ice thickness.  There is an evolving calving front using the eigencalving parameterization (Levermann et al., 2012).  See `README.m` and `run_prog.sh` in the subdirectory.

* `fracture/`: Use PISM's SSA stress balance to do a prognostic run, adding feedback of evolving fracture density on flow, but with a prescribed calving front (`-calving ocean_kill`).  See `README.m` and `run_frac.sh` in the subdirectory.

For each of the three examples an output `.nc` file is produced.  View it with `ncview` or other NetCDF viewer.  If python tools `numpy`, `matplotlib`, and `netcdf4-python` are present, do this

    $ ../plot.py result.nc

from the subdirectory and view the output `.png` images.

Notes
====

* Be aware of periodic boundary conditions at boundaries of the computational domain, and that several fields have jumps at domain boundaries.
* Evolving runs need stress boundary conditions at the ice shelf front (`-cfbc`) plus removal of icebergs if they appear.  Thus you can do runs like this with `-pik` or `-cfbc -part_grid -kill_icebergs`, but simpler option combinations `-cfbc` alone or `-cfbc -part_grid` or `-cfbc -part_grid` or `-cfbc -kill_icebergs` may fail.
