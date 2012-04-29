Ross flow model example
=================

This example demonstrates regional modeling of ice shelves.  We do only a part of a time-evolving model for the changing geometry of the ice shelf.  Namely we only do the equivalent of one-time step, a "diagnostic run" to compute flow velocity for the ice shelf itself.  Therefore this example is an update of the EISMINT-Ross intercomparison experiment (below), but using more modern data.  In particular, we use 5km data for the entire Antarctic ice sheet, therefore this example can be followed to build models of any Antarctic ice shelf.

These publications show the results of PISM applications to Antarctic ice shelves:
* T. Albrecht, A. Levermann (2012).  *Fracture field for large-scale ice dynamics*.  Journal of Glaciology 58 (207), 165--176.
* A. Levermann, T. Albrecht, R. Winkelmann, M. A. Martin, M. Haseloff, I. Joughin (2012) *Kinematic first-order calving law implies potential for abrupt ice-shelf retreat*.  The Cryosphere 6, 273--286. (http://www.the-cryosphere.net/6/273/2012/tc-6-273-2012.html)
* M. A. Martin, R. Winkelmann, M. Haseloff, T. Albrecht, E. Bueler, C. Khroulev, A. Levermann (2011).  *The Potsdam Parallel Ice Sheet Model (PISM-PIK) – Part 2: Dynamic equilibrium simulation of the Antarctic ice sheet*. The Cryosphere 5, 727--740. (http://www.the-cryosphere.net/5/727/2011/tc-5-727-2011.pdf)

The early EISMINT-Ross intercomparison is documented by
* D. R. MacAyeal, V. Rommelaere, Ph. Huybrechts, C.L. Hulbe, J. Determann, C. Ritz (1996). *An ice-shelf model test based on the Ross ice shelf*.  Ann. Glaciol. 23, 46--51. (http://homepages.vub.ac.be/~phuybrec/pdf/MacAyeal.Ann.Glac.23.pdf)
and it appears as an example in the PISM User's Manual in versions 0.4 and earlier.

This directory contains two scripts:
* `preprocess.py`: downloads and preprocessing input data, and
* `run.sh`: use PISM's SSA stress balance to do a "diagnostic" run which computes flow velocity from geometry, ice hardness (a function of temperature), and observed velocities at the grounding line as boundary conditions.

Input datasets
==========

We use the two public datasets:
* *An improved Antarctic dataset for high resolution numerical ice sheet models (ALBMAP v1)*  by A. M. Le Brocq, A. J. Payne, and A. Vieli (http://doi.pangaea.de/10.1594/PANGAEA.734145)
* *MEaSUREs InSAR-Based Antarctica Velocity Map*  by Rignot, E., J. Mouginot, and B. Scheuchl. 2011, (http://nsidc.org/data/nsidc-0484.html)

Fields used
========

* ice thickness (equivalently: the upper/lower surface elevation pair),
* ice surface temperature,
* ice surface mass balance,
* bedrock topography, and
* u and v components of the surface ice velocity for grounded ice as boundary conditions for the floating ice shelf.

Notes
====

* periodic boundary conditions at boundaries of the modeling domain
* several fields have jumps at domain boundaries
* need boundary conditions at the ice shelf front (`-cfbc`)







