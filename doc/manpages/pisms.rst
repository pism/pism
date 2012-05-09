.. The manual page name has to go first, as a top-level header.

=====
pisms
=====

.. The first sub-section header should contain the one-line description

------------------------------------------
PISM's simplified geometry mode executable
------------------------------------------

.. The following are needed to specify the manual page section, group, etc. This seems to be the only way.

:Author: ckhroulev@alaska.edu
:Date:   2011-5-12
:Copyright: Copyright (C) 2011, 2012 Constantine Khroulev
:Version: 0.1
:Manual section: 1
:Manual group: science

SYNOPSIS
========

|  pisms -eisII *experiment name* ...
|  pisms -eisII *experiment name* -i *file.nc* ...
|  pisms -pst ...
|  pisms -mismip ...

DESCRIPTION
===========

**pisms** implements EISMINT II experiments [EISMINT]_, plastic till experiements described in [BBSSA]_ and the [MISMIP]_ setup. Please see the *User's Manual* for details.

.. [EISMINT] **A. Payne**, 1997. EISMINT: Ice sheet model intercomparison exercise phase two. Proposed simplified geometry experiments. http://homepages.vub.ac.be/~phuybrec/eismint/thermo-descr.pdf.

.. [BBSSA] **E. Bueler and J. Brown**, 2009. Shallow shelf approximation as a "sliding law" in a thermodynamically coupled ice sheet model. J. Geophys. Res. 114. F03008, doi:10.1029/2008JF001179.

.. [MISMIP] **C. Schoof, R. Hindmarsh, and F. Pattyn**, 2008. Marine Ice Sheet Model Intercomparison Project. http://homepages.ulb.ac.be/~fpattyn/mismip/.

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

Another important PISM executables include **pismr**, **pismv** and **pross**. 
