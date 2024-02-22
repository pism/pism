.. The manual page name has to go first, as a top-level header.

=====
pismv
=====

.. The first sub-section header should contain the one-line description

---------------------------------------------------
Parallel Ice Sheet Model (PISM): verification tests
---------------------------------------------------

.. The following are needed to specify the manual page section, group, etc. This seems to be the only way.

:Author: ckhroulev@alaska.edu
:Date:   2024-2-21
:Copyright: Copyright (C) 2024 Constantine Khrulev
:Version: 0.1
:Manual section: 1
:Manual group: science

SYNOPSIS
========

|  pismv -test *test_name* -Mx ... -My ... -Mz ... -y ...

DESCRIPTION
===========

**pismv** runs PISM's verification tests. Please see the *PISM User's Manual* for details.

OPTIONS
=======

-i          input file
-y          run length, in model years
-o          output file name
-Mx         number of grid points in the *x* direction
-My         number of grid points in the *y* direction
-Mz         number of grid points in the vertical direction (above the base of the ice)
-Mbz        number of grid points in the vertical direction (in the bed below the base of the ice)
-help       prints PISM and PETSc command-line option help; use with **grep**
-verbose    selects stdout verbosity level, 1 -- no output, 2 -- normal, 3 -- more debugging info, ...

SEE ALSO
========

- The *User's Manual* and other documentation online at https://www.pism.io/docs/manual

- **pismv** to run some verification tests
