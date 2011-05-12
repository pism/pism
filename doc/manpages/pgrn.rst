.. The manual page name has to go first, as a top-level header.

====
pgrn
====

.. The first sub-section header should contain the one-line description

------------------------------
PISM's EISMINT-Greenland setup
------------------------------

.. The following are needed to specify the manual page section, group, etc. This seems to be the only way.

:Author: ckhroulev@alaska.edu
:Date:   2011-5-12
:Copyright: Copyright (C) 2011 Constantine Khroulev
:Version: 0.1
:Manual section: 1
:Manual group: science

SYNOPSIS
========

|  pgrn -i *file.nc* ...
|  pgrn -boot_file *file.nc* ...

DESCRIPTION
===========

**pgrn** implements PISM's [EISGRN]_ setup. This executable is equivalent to **pismr** in every respect except one: it adds an implementation of the EISMINT-Greenland near-air surface temperature parameterization. Please see the *User's Manual* for details.

.. [EISGRN] **C. Ritz**, 1997. EISMINT Intercomparison Experiment: Comparison of existing Greenland models. http://homepages.vub.ac.be/~phuybrec/eismint/greenland.html.

OPTIONS
=======

-i          input file
-y          run length, in model years
-o          output file name
-help       prints PISM and PETSc command-line option help; use with **grep**
-verbose    selects stdout verbosity level, 1 -- no output, 2 -- normal, 3 -- more debugging info, ...

SEE ALSO
========

- The *User's Manual* and other documentation online at http://www.pism-docs.org/
- PISM's *Cheat Sheet* for a quick overview of command-line options.
- PISM's *Source Code Browser* for technical documentation.

Another important PISM executables include **pismr**, **pisms**, **pismv** and **pross**. 
