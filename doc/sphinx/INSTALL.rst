Installing PISM
+++++++++++++++

.. contents:: Table of Contents

The fastest path to a fully functional PISM installation is to use a Linux system with a Debian-based package system (e.g. Ubuntu_): Start by following subsection `Installing prerequisites from Debian packages`_, then `Building PETSc`_, then `Building PISM`_ to install PISM itself.

Libraries and programs needed by PISM
=====================================

This table lists required dependencies for PISM alphabetically.

.. csv-table::
   :header: Required Library, Comment

   FFTW_,       version >= 3.1
   GSL_,        version >= 1.15
   MPI_,        any recent version
   NetCDF_,     version >= 4.1
   PETSc_ [1]_, version >= 3.5
   UDUNITS_,    any recent version

Before installing these "by hand", check sections `Installing prerequisites from Debian packages`_ and `Installing prerequisites using MacPorts on Mac OS X`_ for specific how-to. In particular, if multiple MPI implementations (e.g. MPICH and Open-MPI) are installed then PETSc can under some situations "get confused" and throw MPI-related errors. Even package systems have been known to allow this confusion.

Optional libraries listed below are needed for certain PISM features, namely cell-area correction and parallel I/O. These libraries are recommended, but not strictly required:

.. csv-table::
   :header: Recommended Library, Comment

   PROJ.4_,  Used to compute cell areas and cell bounds
   PnetCDF_, Can be used for parallel I/O


Python_ is needed both in the PETSc installation process and in scripts related to PISM pre- and post-processing, while Git_ is usually needed to download the PISM code. Both should be included in any Linux/Unix distribution.

The following Python packages are needed to do all the examples in the *User’s Manual* (which run Python scripts):

.. csv-table:: Python packages
   :header: Library, Comment

   NumPy_,          used in *most* scripts
   matplotlib_,     used in some scripts
   netcdf4-python_, used in *most* scripts

Installation Cookbook
=====================

Installing PISM requires getting its `prerequisites <Libraries and programs needed by PISM_>`_ and then `building PISM itself <Building PISM_>`_. Sections below cover some cases; please `e-mail us <pism-email_>`_ if you need help.

Installing prerequisites from Debian packages
---------------------------------------------

You should be able to use your package manager to get the prerequisites
for PISM. Install the following packages using ``apt-get`` or
``synaptic`` or similar. All of these are recommended as they satisfy
requirements for building or running PISM.

.. csv-table:: Debian packages
   :header: Package name, Comments

   ``cmake``,            required to configure PISM
   ``libfftw3-dev``,     required by PISM
   ``g++``,              required to build PISM
   ``libgsl0-dev``,      required by PISM
   ``netcdf-bin``,       required: ``ncgen`` is used during the build process
   ``libnetcdf-dev``,    required by PISM
   ``libudunits2-dev``,  required by PISM
   ``cdo``,              used in some pre-processing scripts
   ``cmake-curses-gui``, a text-based easy interface for CMake
   ``git``,              used to get PISM source code
   ``nco``,              used in many pre-processing scripts
   ``ncview``,           view fields in NetCDF files
   ``libproj-dev``,      used to compute ice area and volume
   ``python-dev``,       (helps with scripts…perhaps not essential)
   ``python-pyproj``,    used in some pre-processing scripts
   ``python-netcdf4``,   used in most post-processing scripts
   ``libx11-dev``,       X windows is useful to get graphics through PETSc
   ``libblas-dev``,      BLAS is required by PETSc
   ``liblapack-dev``,    LAPACK is required by PETSc
   ``openmpi-bin``,      MPI is required to run PISM in parallel
   ``libopenmpi-dev``,   MPI is required to run PISM in parallel

You may be able to install these by running

.. code::

   sudo apt-get install cmake libfftw3-dev g++ libgsl0-dev netcdf-bin \
                        libnetcdf-dev libudunits2-dev cdo cmake-curses-gui \
                        git nco ncview libproj-dev python-dev python-pyproj \
                        python-netcdf4 libx11-dev libblas-dev liblapack-dev \
                        openmpi-bin libopenmpi-dev

(You may need to change this command to match your package system.)

Once done, see `Building PETSc`_ to install PETSc from source and then `Building PISM`_ for building PISM itself.

Installing prerequisites using MacPorts on Mac OS X
---------------------------------------------------

Follow these steps to install PISM's prerequisites on the Mac OS X
operating system.

#. As PISM is distributed as source code only, you will need software developer’s tools, XCode_ and the *X window system server*, XQuartz_.

#. The use of MacPorts_ (or Fink_, or Homebrew_) is recommended, as it significantly simplifies installing many open-source libraries. These instructions assume that you use MacPorts_. Download a package from the MacPorts_, install, and set the environment:

   .. code:: bash

      export PATH=/opt/local/bin:/opt/local/sbin:$PATH

   for MacPorts.

#. It may not be necessary to install Python, as it is bundled with the operating system. Some PISM scripts use SciPy; it can be installed using MacPorts or by downloading the `Enthought Python Distribution <Enthought_>`_.

#. If you are using MacPorts, do

   .. code:: bash

      sudo port install git cmake fftw-3 gsl mpich-default netcdf udunits2 libproj4 ncview

#. At this point, all the PISM prerequisites except PETSc are installed. Follow instructions in `Building PETSc`_ to install it.

#. Now you can build PISM as described in section `Building PISM`_.

Installing prerequisites from sources
-------------------------------------

From now on, this manual assumes the use of the Bash_ shell.

#. You will need Python_ and Git_ installed. To use the (recommended) graphical output of PISM you will need an `X Window server <X_>`_.

#. Generally the "header files" for its prerequisite libraries are required for building PISM. (This means that the "developer’s versions" of the libraries are needed if the libraries are downloaded from package repositories like Debian's; see the `prerequisite list <Libraries and programs needed by PISM_>`_.)

#. PISM uses `NetCDF <NetCDF_>`_ as an input and output file format. If it is not already present, install it using the instructions at the web-page or using a package management system.

#. PISM uses the `GNU Scientific Library <GSL_>`_ for certain numerical calculations and special functions. If it is not already present, install it using the instructions at the web-page or using a package management system.

#. PISM uses the `FFTW library <FFTW_>`_ for the deformation of the solid earth (bed) under ice loads. Install FFTW version 3.1 or later, or check that it is installed already.

#. You will need a version of `MPI <MPI_>`_. Your system may have an existing MPI installation, in which case it should probably be used when `building PETSc <Building PETSc_>`_. The goal is to have the PETSc installation use the same version of MPI which is called by the ``mpiexec`` or ``mpirun`` executable.

   If you had to install an MPI library "by hand" you will want to add
   the MPI ``bin`` directory to your path so that you can run parallel
   programs using the ``mpiexec`` or ``mpirun`` command. For example,
   you can add it with the statement

   .. code:: bash

      export PATH=/home/user/mympi/bin:$PATH

   (for Bash shell).

   Such a statement can, of course, appear in your ``.bashrc`` (or
   ``.profile``) file so that there is no need to retype it each time
   you use MPI.

#. PISM uses UDUNITS_ to convert units of physical quantities read from input files and written to output files. Follow instructions on its website to install.

Building PETSc
--------------

PISM is built on top of PETSc_, which is actively developed and an up-to-date PETSc distribution is unlikely to be available in package repositories. Download the PETSc source by grabbing the current gzipped tarball at:

http://www.mcs.anl.gov/petsc/download/index.html

(See `PISM's prerequisite list <Libraries and programs needed by PISM_>`_ for the minimum supported PETSc version.) The "lite" form of the tarball is fine if you are willing to depend on an Internet connection for accessing PETSc documentation.

You should configure and build PETSc as described on the PETSc installation page, but it might be best to read the following comments on the PETSc configure and build process first:

#. Untar in your preferred location and enter the new PETSc directory. Note PETSc should *not* be configured using root privileges. When you run the configure script the following options are recommended; note PISM uses shared libraries by default:

   .. code:: bash

      export PETSC_DIR=$PWD
      export PETSC_ARCH=opt
      ./config/configure.py --with-shared-libraries --with-debugging=0 --with-fc=0

   You need to define the environment variables ``PETSC_DIR`` and ``PETSC_ARCH`` [6]_ -- one way is shown here-- *before* running the configuration script. Turning off the inclusion of debugging code and symbols can give a significant speed improvement, but some kinds of development will benefit from setting ``--with-debugging=1``. Using shared libraries may be unwise on certain clusters; check with your system administrator. PISM does not use PETSc's Fortran API, so the Fortran compiler is disabled by ``--with-fc=0``.

#. It is sometimes convenient to have PETSc grab a local copy of BLAS and LAPACK rather than using the system-wide version. So one may add "``--download-f2cblaslapack=1``" to the other configure options.

#. If there is an existing MPI installation, for example at ``/home/user/mympi/``, one can point PETSc to it by adding the option "``--with-mpi-dir=/home/user/mympi/``". The path used in this option must have MPI executables ``mpicxx`` and ``mpicc``, and either ``mpiexec`` or ``mpirun``, in sub-directory ``bin/`` and MPI library files in sub-directory ``lib/``. Alternatively, use MPI's compiler wrappers to specify an MPI library when installing PETSc, for example:

   .. code:: shell

      CC=mpicc CXX=mpicxx ./config/configure.py --with-shared-libraries --with-debugging=0 --with-fc=0

   If you get messages suggesting that PETSc cannot configure using your existing MPI, you might want to try adding the ``--download-mpich=1`` (or ``--download-openmpi=1``) option to PETSc’s configure command.

#. Configuration of PETSc for a batch system requires special procedures described at the PETSc documentation site. One starts with a configure option ``--with-batch=1``. See the "Installing on machine requiring cross compiler or a job scheduler" section of the `PETSc installation page <PETSc-installation_>`_.

#. Configuring PETSc may take a moment even when everything goes smoothly. A value for the environment variable ``PETSC_ARCH`` will be reported at the end of the configure process; take note of this value. One may always reconfigure with additional ``PETSC_ARCH`` as needed.

#. After ``configure.py`` finishes, you will need to ``make all test`` in the PETSc directory and watch the result. If the X Windows system is functional some example viewers will appear; as noted you will need the X header files for this to work.

Building PISM
-------------

At this point you have configured the environment which PISM needs.

To make sure that the key PETSc and MPI prerequisites work properly together, so that you can run PISM in parallel, you might want to make sure that the correct ``mpiexec`` can be found, by setting your ``PATH``. For instance, if you used the option ``--download-mpich=1`` in the PETSc configure, the MPI ``bin`` directory will have a path like ``$PETSC_DIR/$PETSC_ARCH/bin``. Thus the following lines might appear in your ``.bashrc`` or ``.profile``, if not there already:

.. code:: bash

   export PETSC_DIR=/home/user/petsc-3.7.0/
   export PETSC_ARCH=opt
   export PATH=$PETSC_DIR/$PETSC_ARCH/bin/:$PATH

From now on we will assume that the ``PETSC_ARCH`` and ``PETSC_DIR`` variables are set.

You are ready to build PISM itself, which is a much quicker procedure, as follows:

#. Get the latest source for PISM using the Git_ version control system:

   Check `PISM's website <PISM_>`_ for the latest version of PISM.

   .. _git-clone:

   Run

   .. code:: bash

      git clone git://github.com/pism/pism.git pism-stable

   A directory called "``pism-stable``" will be created. Note that in the future when you enter that directory, ``git pull`` will update to the latest revision of PISM. [2]_

#. Build PISM:[3]_

   .. code:: bash

      mkdir -p pism-stable/build
      cd pism-stable/build
      PISM_INSTALL_PREFIX=~/pism cmake ..
      make install

   Here ``pism-stable`` is the directory containing PISM source code while ``~/pism`` is the directory PISM will be installed into. All the temporary files created during the build process will be in ``pism-stable/build`` created above.

   You might need to add ``CC`` and ``CXX`` to the ``cmake`` command:

   .. code:: bash

      PISM_INSTALL_PREFIX=~/pism CC=mpicc CXX=mpicxx cmake ..

   Whether this is necessary or not depends on your MPI setup.

   Commands above will configure PISM to be installed in ``~/pism/bin`` and ``~/pism/lib/`` then compile and install all its executables and scripts.

   If your operating system does not support shared libraries [4]_, then set ``Pism_LINK_STATICALLY`` to "ON". This can be done by either running

   .. code:: bash

      cmake -DPism_LINK_STATICALLY=ON ..

   or by using ``ccmake`` [5]_ run

   .. code:: bash

      ccmake ..

   and then change ``Pism_LINK_STATICALLY`` (and then press ``c`` to "configure" and ``g`` to "generate Makefiles"). Then run ``make install``.

   Object files created during the build process (located in the ``build`` sub-directory) are not automatically deleted after installing PISM, so run "``make clean``" if space is an issue. You can also delete the build directory altogether if you are not planning on re-compiling PISM.

   Note that when using Intel's compiler high optimization settings such as ``-O3``, ``-fp-model precise`` may be needed to get reproducible model results. Set it using ``ccmake`` or by setting ``CFLAGS`` and ``CXXFLAGS`` environment variables when building PISM's prerequisites and PISM itself.

   .. code:: bash

      export CFLAGS="-fp-model precise"
      export CXXFLAGS="-fp-model precise"
      cmake [other options] ..

#. PISM executables can be run most easily by adding the ``bin/`` sub-directory in your selected install path (``~/pism/bin`` in the example above) to your ``PATH``. For instance, this command can be done in the Bash_ shell or in your ``.bashrc`` file:

   .. code:: bash

      export PATH=~/pism/bin:$PATH

#. Now see section `Quick tests of the installation`_ or the *Getting Started* section of the *User’s Manual* to continue.

Common build problems and solutions
-----------------------------------

We recommend using ``ccmake``, the text-based CMake interface to adjust PISM’s build parameters. One can also set CMake cache variables using the ``-D`` command-line option (``cmake -Dvariable=value``) or by editing ``CMakeCache.txt`` in the build directory.

Here are some issues we know about.

-  Sometimes, if a system has more than one MPI installation CMake finds the wrong one. To tell it which one to use, set ``MPI_LIBRARY`` and related variables by using ``ccmake``. You can also set environment variables ``CC`` and ``CXX`` to point to MPI wrappers:

   .. code:: bash

      CC=mpicc CXX=mpicxx cmake path/to/pism-source

   It is also possible to guide CMake’s configuration mechanism by setting ``MPI_COMPILER`` to the compiler (such as ``mpicc``) corresponding to the MPI installation you want to use, setting ``MPI_LIBRARY`` to ``MPI_LIBRARY-NOTFOUND`` and re-running CMake.

-  If you are compiling PISM on a system using a cross-compiler, you will need to disable CMake’s tests trying to determine if PETSc is installed properly. To do this, set ``PETSC_EXECUTABLE_RUNS`` to "yes".

   To tell CMake where to look for libraries for the target system, see `CMake cross compiling <CMake-cross-compiling_>`_ and the paragraph about ``CMAKE_FIND_ROOT_PATH`` in particular.

-  Note that the PISM build system uses ``ncgen`` from the NetCDF package to generate the configuration file ``pism_config.nc``. This means that a working NetCDF installation is required on both the "host" and the "target" systems when cross-compiling PISM.

-  Some systems support static libraries only. To build PISM statically and tell CMake not to try to link to shared libraries, set ``Pism_LINK_STATICALLY`` to ``ON`` using ``ccmake``.

-  You can set ``Pism_LOOK_FOR_LIBRARIES`` to "``OFF``" to disable all heuristics and set compiler flags by hand. See `HPC builds <HPC-builds_>`_ for examples.

Quick tests of the installation
===============================

Once you’re done with the installation, a few tests can confirm that PISM is functioning correctly.

#. Try a MPI four process verification run:

   .. code:: bash

      mpiexec -n 4 pismv -test G -y 200

   If you see some output and a final ``Writing model state`` ``to file ’unnamed.nc’`` then PISM completed successfully. At the end of this run you get measurements of the difference between the numerical result and the exact solution. See the *User’s Manual* for more on PISM verification.

   The above "``-n 4``" run should work even if there is only one actual processor (core) on your machine. (In that case MPI will just run multiple processes on the one processor.) This run will also produce a NetCDF output file ``unnamed.nc``, which can be read and viewed by NetCDF tools.

#. Try an EISMINT II run using the PETSc viewers (under the X window system):

   .. code::

      pisms -y 5000 -view thk,temppabase,velsurf_mag

   When using such viewers and ``mpiexec`` the additional final option ``-display :0`` is sometimes required to enable MPI to use X, like this:

   .. code::

       mpiexec -n 2 pisms -y 5000 -view thk,temppabase,velsurf_mag -display :0

   Also ``-drawpause 0.1`` or similar may be needed if the figures are refreshing too fast.

#. Run a basic suite of software tests. To do this, make sure that NCO_ and Python packages NumPy_ and netcdf4-python_ are installed. Also, the CMake flag ``Pism_BUILD_EXTRA_EXECS`` should be ``ON``. Then run:

   .. code:: bash

      make       # do this if you changed something with CMake
      make test

   in the build directory. The message at the bottom should say "``100% tests passed, 0 tests failed out of XX``" or similar. Feel free to `send us <pism-email_>`_ the output of ``make test``. if any failed tests cannot be resolved.

Next steps
==========

Start with the *User’s Manual*, which has a "Getting started" section. A copy is on-line at the `PISM homepage <PISM_>`_, along with a `source code <pism-code-browser_>`_ (HTML). Completely up-to-date documentation can be built from LaTeX source in the ``doc/`` sub-directory, as described in the next section.

A final reminder with respect to installation: Let’s assume you have checked out a copy of PISM using Git, `as described above <git-clone_>`_. You can then update your copy of PISM to the latest version by running ``git pull`` in the PISM directory and ``make install`` in your build directory.

Rebuilding PISM documentation
=============================

You might want to rebuild the documentation from source, as PISM and its
documentation evolve together. These tools are required:

.. csv-table::
   :header: Tool, Comment

   LaTeX_,    needed for rebuilding any of the documentation
   doxygen_,  required to rebuild the *Browser* from source
   graphviz_, required to rebuild the *Browser* from source

To rebuild PISM documentation, change to the PISM build directory and do

.. csv-table::
   :header: Command, Comment

   ``make pism_manual``,  "to build the *User’s Manual*, ``pism_manual.pdf``"
   ``make pism_forcing``, "to build the *PISM’s Climate Forcing Components* document, ``pism_forcing.pdf``"
   ``make browser``,       to build the *PISM Source Code Browser*.

To build documentation on a system without PISM’s prerequisite libraries
(such as MPI and PETSc), assuming that PISM sources are in
``~/pism-stable``, do the following:

.. code:: shell

   cd ~/pism-stable
   mkdir doc-build # create a build directory
   cd doc-build
   cmake ../doc

then commands "``make pism_manual``", "``make pism_forcing``" and others (see above) will work as expected.

Building documentation for PISM’s Python bindings and inversion tools
---------------------------------------------------------------------

The documentation for PISM’s Python bindings uses the documentation-generation tool Sphinx_. The bindings make scripting and interactive PISM possible, but many PISM users will not need them. Installing them is required to use PISM for inversion of surface velocities for basal shear stress and ice hardness. Building their documentation is strongly-recommended before use.

Sphinx_ can be installed using ``apt-get`` or MacPorts_; see the website for more details. For example, do

.. code:: shell

   sudo apt-get install sphinx-common

The bindings documentation also requires the Sphinx extension called ``sphinxcontrib.bibtex``, which may come with some Sphinx packages (but not with Debian packages at this time). Without it you will see this error when you try to build the bindings documentation:

.. code::

   Extension error:
   Could not import extension sphinxcontrib.bibtex (exception: No module named bibtex)

To install it see http://sphinxcontrib-bibtex.readthedocs.io/en/latest/.

Note that if you install Sphinx using MacPorts_, you will install a version that depends on your Python version, and its executables will have names that depend on the Python version, e.g. ``sphinx-build-2.7`` rather than ``sphinx-build`` for Python 2.7. You will want to set up aliases so that the standard names work as well. To do this,

.. code::

    sudo port select sphinx py27-sphinx

(replacing ``py27-sphinx`` with ``py26-sphinx`` for Python 2.6, etc.) If you opt not to do this, you can tell CMake the name of your Sphinx executable using

.. code::

   cmake -DSPHINX_EXECUTABLE=sphinx-build-2.7 ...

for example.

Now you can build the documentation. In the PISM build directory, do

.. code::

    make pismpython_docs

If you get an error like

.. code::

   make: *** No rule to make target `pismpython_docs'.  Stop.

then re-run ``cmake ..`` or ``ccmake ..``, making sure that Sphinx is installed (see above); the ``pismpython_docs`` make target will then be present.

The main page for the documentation is then in ``doc/pismpython/html/index.html`` inside your build directory. The documentation build can take some time while it builds a large number of small images from LaTeX formulas.

.. rubric:: Footnotes

.. [1] "PETSc" is pronounced "pet-see".

.. [2] Of course, after ``git pull`` you will ``make -C build install`` to recompile and re-install PISM.

.. [3] Please report any problems you meet at these build stages by `sending us <pism-email_>`_ the output.

.. [4] This might be necessary if you’re building on a Cray XT5 or a Sun Opteron Cluster, for example.

.. [5] Install the ``cmake-curses-gui`` package to get ``ccmake`` on Ubuntu_.

.. [6] The ``PETSC_ARCH`` variable is just a string you can use to choose different PETSc configurations and does not have any other significance.

.. _Bash: http://www.gnu.org/software/bash/
.. _CMake-cross-compiling: https://cmake.org/Wiki/CMake_Cross_Compiling
.. _doxygen: http://www.stack.nl/~dimitri/doxygen/
.. _Enthought: https://www.enthought.com/
.. _FFTW: http://www.fftw.org/
.. _Fink: http://www.finkproject.org/
.. _GSL: http://www.gnu.org/software/gsl/
.. _Git: https://git-scm.com/
.. _graphviz: http://www.graphviz.org/
.. _Homebrew: https://brew.sh/
.. _LaTeX: http://www.latex-project.org/
.. _MPI: http://www.mcs.anl.gov/research/projects/mpi/
.. _MacPorts: https://www.macports.org/
.. _NCO: http://nco.sourceforge.net/
.. _NetCDF: http://www.unidata.ucar.edu/software/netcdf/
.. _NumPy: http://www.numpy.org/
.. _PETSc: http://www.mcs.anl.gov/petsc/
.. _PETSc-installation: http://www.mcs.anl.gov/petsc/documentation/installation.html
.. _PROJ.4: http://proj4.org/
.. _PnetCDF: http://trac.mcs.anl.gov/projects/parallel-netcdf
.. _Python: https://www.python.org
.. _Sphinx: http://www.sphinx-doc.org/en/stable/install.html
.. _UDUNITS: http://www.unidata.ucar.edu/software/udunits/
.. _Ubuntu: https://www.ubuntu.com/desktop
.. _X: https://www.x.org/wiki/
.. _XCode: https://developer.apple.com/xcode/
.. _XQuartz: https://www.xquartz.org/
.. _matplotlib: http://matplotlib.org/
.. _netcdf4-python: https://pypi.python.org/pypi/netCDF4

.. _HPC-builds: https://github.com/pism/hpc-builds
.. _pism-email: mailto:uaf-pism@alaska.edu
.. _PISM: http://www.pism-docs.org/wiki/doku.php
.. _pism-code-browser: http://www.pism-docs.org/doxy/html/index.html

..
   Local Variables:
   eval: (visual-line-mode nil)
   fill-column: 1000
   End:
