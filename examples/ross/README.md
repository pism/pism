Ross flow model example
=================

This example demonstrates regional modeling of ice shelves.  We do only one-time step, a "diagnostic run" to compute the flow velocity of the floating ice.  This example is an update of the EISMINT-Ross intercomparison experiment (below), using more recent data, a 5km data set for the entire Antarctic ice sheet.  This example can be modified to build diagnostic models of any Antarctic ice shelf.

These publications show the results of PISM applications to Antarctic ice shelves:
* T. Albrecht, A. Levermann (2012).  *Fracture field for large-scale ice dynamics*.  Journal of Glaciology 58 (207), 165--176.
* A. Levermann, T. Albrecht, R. Winkelmann, M. A. Martin, M. Haseloff, I. Joughin (2012) *Kinematic first-order calving law implies potential for abrupt ice-shelf retreat*.  The Cryosphere 6, 273--286. (http://www.the-cryosphere.net/6/273/2012/tc-6-273-2012.html)
* M. A. Martin, R. Winkelmann, M. Haseloff, T. Albrecht, E. Bueler, C. Khroulev, A. Levermann (2011).  *The Potsdam Parallel Ice Sheet Model (PISM-PIK) – Part 2: Dynamic equilibrium simulation of the Antarctic ice sheet*. The Cryosphere 5, 727--740. (http://www.the-cryosphere.net/5/727/2011/tc-5-727-2011.pdf)

The early EISMINT-Ross intercomparison, which appears as an example in the PISM User's Manual in versions 0.4 and earlier, is documented by
* D. R. MacAyeal, V. Rommelaere, Ph. Huybrechts, C.L. Hulbe, J. Determann, C. Ritz (1996). *An ice-shelf model test based on the Ross ice shelf*.  Ann. Glaciol. 23, 46--51. (http://homepages.vub.ac.be/~phuybrec/pdf/MacAyeal.Ann.Glac.23.pdf)

Input datasets
==========

We use the two public datasets:
* *An improved Antarctic dataset for high resolution numerical ice sheet models (ALBMAP v1)*  by A. M. Le Brocq, A. J. Payne, and A. Vieli (http://doi.pangaea.de/10.1594/PANGAEA.734145)
* *MEaSUREs InSAR-Based Antarctica Velocity Map*  by Rignot, E., J. Mouginot, and B. Scheuchl. 2011, (http://nsidc.org/data/nsidc-0484.html)

From these data sets, the following fields are used:
* ice thickness (equivalently: the upper/lower surface elevation pair),
* ice surface temperature,
* ice surface mass balance,
* bedrock topography, and
* u and v components of the surface ice velocity for grounded ice as boundary conditions for the floating ice shelf.


Running the example
=================

This directory contains two scripts:
* `preprocess.py`: downloads and preprocessing input data, and
* `run.sh`: use PISM's SSA stress balance to do a "diagnostic" run which computes flow velocity from geometry, ice hardness (a function of temperature), and observed velocities at the grounding line as boundary conditions.

Do:

    $ ./preprocess.py
    $ ./run.sh

Notes
====

* periodic boundary conditions at boundaries of the modeling domain
* several fields have jumps at domain boundaries
* clearly needs boundary conditions at the ice shelf front (`-cfbc`); runs with `-pik` or `-cfbc -part_grid -kill_icebergs` but not with `-cfbc` alone or `-cfbc -part_grid` or `-cfbc -part_grid` or `-cfbc -kill_icebergs`


