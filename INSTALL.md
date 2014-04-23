### Welcome to PISM!

All information about PISM can be found at

[`www.pism-docs.org`](http://www.pism-docs.org)

To download the latest stable release use git:

    git clone git://github.com/pism/pism.git pism0.6

For complete installation instructions get the PISM Installation Manual (PDF)
from http://www.pism-docs.org.  Major installation requirements:

- You must have [MPI](http://www.mcs.anl.gov/mpi/), [GSL](http://www.gnu.org/software/gsl/),
  [NetCDF](http://www.unidata.ucar.edu/software/netcdf/), [FFTW-3](http://www.fftw.org),
  and [UDUNITS-2](http://www.unidata.ucar.edu/software/udunits/).

- [PETSc](http://www-unix.mcs.anl.gov/petsc/) must be installed and functioning
  before building PISM.  Environment variables `PETSC_DIR` and `PETSC_ARCH`
  must be set.  PETSc 3.3 or later is required.

To build PISM locally (within the PISM source tree), do

    cd pism0.6
    mkdir build
    cd build
    ccmake ..     # hit 'c' to see initial config, 'c' again to configure,
                  # and then 'g' to generate once settings look right
    make install

Check where PISM executables can be found by:

    which pismr

Run a simplified geometry experiment on one MPI process:

    pisms

Everything is working if this run goes for 1000 model years and ends with
"... done with run".  Again but with 2 MPI processes:

    mpiexec -n 2 pisms

A more thorough test of the build requires tools (esp. Python, NumPy, and NCO):

    make test

Please see the PISM User's Manual, also at [www.pism-docs.org](http://www.pism-docs.org).
It explains actual ice sheet modeling with PISM.

To generate PDF documentation:

    cd 'your build directory'
    make manual              # builds manual.pdf
    make installation        # builds installation.pdf
    make forcing             # builds forcing.pdf

To generate Doxygen source code documentation:

    make browser             # doxygen docs; view doc/browser/html/index.html when done

To generate PDF and Doxygen documentation on a system without PETSc and other prerequisites, do

    cd pism0.6
    mkdir build-doc
    cd build-doc
    cmake ../doc

and then run `make` as described above (in this build directory).

See also [the development version page](http://www.pism-docs.org/wiki/doku.php?id=development_version).

For questions about PISM usage, installation, creating new models, etc.,
e-mail <help@pism-docs.org>.

