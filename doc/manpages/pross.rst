.. The manual page name has to go first, as a top-level header.

=====
pross
=====

.. The first sub-section header should contain the one-line description

-------------------------
PISM's EISMINT-Ross setup
-------------------------

.. The following are needed to specify the manual page section, group, etc. This seems to be the only way.

:Author: ckhroulev@alaska.edu
:Date:   2011-5-12
:Copyright: Copyright (C) 2011, 2012 Constantine Khroulev
:Version: 0.1
:Manual section: 1
:Manual group: science

SYNOPSIS
========

|  pross -boot_file *file.nc* -Mx ... -My ... -riggs ...

DESCRIPTION
===========

**pross** implements PISM's [EISROSS]_ setup. This is a special executable that does a *diagnostic* run only. It uses PISM's stress balance code *independently* from the code implementing time stepping, providing climate boundary conditions, etc. Please see the *User's Manual* for details.

.. [EISROSS] **D. R. MacAyeal, V. Rommelaere, P. Huybrechts, C. Hulbe, J. Determann, and C. Ritz**, 1996. An ice-shelf model test based on the Ross ice shelf. Ann. Glaciol. 23, 46–51.

OPTIONS
=======

-boot_file  input file
-o          output file name
-help       prints PISM and PETSc command-line option help; use with **grep**
-verbose    selects stdout verbosity level, 1 -- no output, 2 -- normal, 3 -- more debugging info, ...

SEE ALSO
========

- The *User's Manual* and other documentation online at http://www.pism-docs.org/
- PISM's *Source Code Browser* for technical documentation.

Another important PISM executables include **pismr**, **pisms** and **pismv**.
