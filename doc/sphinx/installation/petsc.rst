.. include:: ../global.txt

.. _sec-install-petsc:

Building PETSc
--------------

PISM is built on top of PETSc_, which is actively developed and an up-to-date PETSc
distribution is unlikely to be available in package repositories. Download the PETSc
source by grabbing the current gzipped tarball at:

    |petsc-download|

(Use version |petsc-min-version| or newer; see :ref:`sec-install-prerequisites` for
details.) The "lite" form of the tarball is fine if you are willing to depend on an
Internet connection for accessing PETSc documentation.

You should configure and build PETSc as described on the PETSc installation page, but it
might be best to read the following comments on the PETSc configure and build process
first:

#. Untar in your preferred location and enter the new PETSc directory. Note PETSc should
   *not* be configured using root privileges. When you run the configure script the
   following options are recommended; note PISM uses shared libraries by default:

   .. code-block:: bash

      export PETSC_DIR=$PWD
      export PETSC_ARCH=opt
      ./config/configure.py --with-shared-libraries \
                            --with-debugging=0 \
                            --with-fc=0

   .. note::

      PETSc's ``configure.py`` requires Python 2.x (Python 3.x is not supported yet).

   You need to define the environment variables ``PETSC_DIR`` and ``PETSC_ARCH`` [#]_ --
   one way is shown here -- *before* running the configuration script. Turning off the
   inclusion of debugging code and symbols can give a significant speed improvement, but
   some kinds of development will benefit from setting ``--with-debugging=1``. Using
   shared libraries may be unwise on certain clusters; check with your system
   administrator. PISM does not use PETSc's Fortran API, so the Fortran compiler is
   disabled by ``--with-fc=0``.

#. It is sometimes convenient to have PETSc grab a local copy of BLAS and LAPACK rather
   than using the system-wide version. So one may add "``--download-f2cblaslapack=1``" to
   the other configure options.

#. If there is an existing MPI installation, we recommend using MPI's compiler wrappers to
   specify an MPI library when installing PETSc, for example:

   .. code-block:: bash

      CC=mpicc CXX=mpicxx ./config/configure.py --with-shared-libraries \
                                                --with-debugging=0 \
                                                --with-fc=0

   If you get messages suggesting that PETSc cannot configure using your existing MPI, you
   might want to try adding the ``--download-mpich=1`` (or ``--download-openmpi=1``)
   option to PETSc’s configure command.

#. Configuration of PETSc for a batch system requires special procedures described at the
   PETSc documentation site. One starts with a configure option ``--with-batch=1``. See
   the "Installing on machine requiring cross compiler or a job scheduler" section of the
   `PETSc installation page <PETSc-installation_>`_.

#. Configuring PETSc may take a moment even when everything goes smoothly. A value for the
   environment variable ``PETSC_ARCH`` will be reported at the end of the configure
   process; take note of this value. One may always reconfigure with additional
   ``PETSC_ARCH`` as needed.

#. After ``configure.py`` finishes, you will need to ``make all test`` in the PETSc
   directory and watch the result. If the X Windows system is functional some example
   viewers will appear; as noted you will need the X header files for this to work.

.. rubric:: Footnotes

.. [#] The ``PETSC_ARCH`` variable is just a string you can use to choose different PETSc
       configurations and does not have any other significance.

   
