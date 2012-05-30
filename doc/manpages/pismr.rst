=====
pismr
=====

--------------------------
PISM's run mode executable
--------------------------
:Author: ckhroulev@alaska.edu
:Date:   2011-5-12
:Copyright: Copyright (C) 2011, 2012 Constantine Khroulev
:Version: 0.1
:Manual section: 1
:Manual group: science

SYNOPSIS
========

|  pismr -i *input_file.nc* ...
|  pismr -boot_file *bootstrapping_file.nc* ...

DESCRIPTION
===========

**pismr** performs an evolution PISM modeling run. Please see the *User's Manual* for details.

OPTIONS
=======

-i          input file
-boot_file  specifies a bootstrapping file
-y          run length, in model years
-o          output file name
-help       prints PISM and PETSc command-line option help; use with **grep**
-verbose    selects stdout verbosity level, 1 -- no output, 2 -- normal, 3 -- more debugging info, ...

SEE ALSO
========

- The *User's Manual* and other documentation online at http://www.pism-docs.org/
- PISM's *Source Code Browser* for technical documentation.

Another important PISM executables include **pisms**, **pismv** and **pross**. 
