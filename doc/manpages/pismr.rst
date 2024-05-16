=====
pismr
=====

------------------------------------------------
Parallel Ice Sheet Model (PISM): prognostic runs
------------------------------------------------
:Author: ckhroulev@alaska.edu
:Date:   2024-2-21
:Copyright: Copyright (C) 2024 Constantine Khrulev
:Version: 0.1
:Manual section: 1
:Manual group: science

SYNOPSIS
========

|  pismr -i *input_file.nc* ...
|  pismr -boot_file *bootstrapping_file.nc* ...

DESCRIPTION
===========

``pismr`` performs an evolution run using PISM, a Parallel Ice Sheet Model. Please see the
*PISM User's Manual* for details.

OPTIONS
=======

-i          input file
-bootstrap  enables bootstrapping mode
-y          run length, in model years
-o          output file name
-help       prints PISM and PETSc command-line option help; use with **grep**
-verbose    selects stdout verbosity level, 1 -- minimal output, 2 -- normal, 3 -- more debugging info, ...
-eisII X    run an EISMINT II experiment X
-test X     run a verification test X

SEE ALSO
========

- The *User's Manual* and other documentation online at https://www.pism.io/docs/manual
