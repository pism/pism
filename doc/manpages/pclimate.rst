.. The manual page name has to go first, as a top-level header.

========
pclimate
========

.. The first sub-section header should contain the one-line description

--------------------------------------------
PISM's executable for testing climate inputs
--------------------------------------------

.. The following are needed to specify the manual page section, group, etc. This seems to be the only way.

:Author: ckhroulev@alaska.edu
:Date:   2011-5-12
:Copyright: Copyright (C) 2011, 2012 Constantine Khroulev
:Version: 0.1
:Manual section: 1
:Manual group: science

SYNOPSIS
========

|  pclimate -i *file.nc* -o *output_file.nc* -times *start_year:dt:end_year* *other options*

DESCRIPTION
===========

**pclimate** uses a PISM output file as its input. It reads model state data
from this file, sets up *surface*, *atmosphere* and *ocean* models and saves
boundary conditions provided by these models at specified times to the output
file.

This allows testing PISM's climate inputs that **do not** depend on
ice-dynamics-driven feedback.

Please see the *User's Manual* for details.

OPTIONS
=======

-i  input file name
-o  output file name
-times  spcifies times at which to save climate data
-surface  selects a surface model
-atmosphere  selects an atmosphere model
-ocean  selects an ocean model
-help  prints PISM and PETSc command-line option help; use with **grep**
-verbose  selects stdout verbosity level, 1 -- no output, 2 -- normal, 3 -- more debugging info, ...

SEE ALSO
========

- The *User's Manual* and other documentation online at http://www.pism-docs.org/
- PISM's *Source Code Browser* for technical documentation.

Another important PISM executables include **pismr**, **pisms** and **pismv**.
