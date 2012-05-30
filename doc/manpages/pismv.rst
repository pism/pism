.. The manual page name has to go first, as a top-level header.

=====
pismv
=====

.. The first sub-section header should contain the one-line description

-----------------------------------
PISM's verification mode executable
-----------------------------------

.. The following are needed to specify the manual page section, group, etc. This seems to be the only way.

:Author: ckhroulev@alaska.edu
:Date:   2011-5-12
:Copyright: Copyright (C) 2011, 2012 Constantine Khroulev
:Version: 0.1
:Manual section: 1
:Manual group: science

SYNOPSIS
========

|  pismv -test *test_name* -Mx ... -My ... -Mz ... -y ...

DESCRIPTION
===========

**pismv** implements PISM's verification tests. Please see the *User's Manual* for details.

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
- PISM's *Source Code Browser* for technical documentation.

Another important PISM executables include **pismr**, **pisms** and **pross**. 
