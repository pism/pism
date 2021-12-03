.. include:: ../global.txt

.. _sec-install-petsc:

Building PETSc
--------------

PISM is built on top of PETSc_, which is actively developed and an up-to-date PETSc
distribution may not be available in package repositories. Download the PETSc
source by grabbing the current gzipped tarball at:

    |petsc-download|

Use version |petsc-min-version| or newer; see :ref:`sec-install-prerequisites` for
details.

You should configure and build PETSc as described on the `PETSc installation page
<PETSc-installation_>`_, but it might be best to read the following comments on the PETSc
configure and build process first.

Untar in your preferred location and enter the new PETSc directory. Note PETSc should
*not* be configured using root privileges. When you run the configure script the following
options are recommended; note PISM uses shared libraries by default:

.. literalinclude:: code/petsc.sh
   :language: bash
   :caption: Building PETSc
   :start-after: manual-begin
   :end-before: manual-end

This will install PETSc and its Python bindings in the directory ``~/local/petsc``. Remove
``--with-petsc4py`` if you don't need Python bindings (e.g. if you are not going to use
PISM's Python bindings).

- You need to define the environment variables ``PETSC_DIR`` and ``PETSC_ARCH``\ [#f1]_
  (one way is shown here) *before* running the configuration script.

- We recommend using MPI's compiler wrappers to specify an MPI library when installing
  PETSc (see ``--with-cc=mpicc``, ``--with-cxx=mpicxx``, and ``--with-fc=mpifort`` above).

- Turning off the inclusion of debugging code and symbols (``--with-debugging=0``) can
  give a significant speed improvement, but some kinds of development will benefit from
  setting ``--with-debugging=1``.

- Using shared libraries may be unwise on certain clusters; check with your system
  administrator.

- The option ``--download-f2cblaslapack=1`` tells PETSc to download BLAS and LAPACK rather
  than using the system-wide version. This tends to work well for PISM, but see section
  3.5.3 of :cite:`petsc-user-ref` for instructions regarding building PETSc with optimized
  BLAS and LAPACK libraries.

- If you get messages suggesting that PETSc cannot configure using your existing MPI, you
  might want to try adding ``--download-mpich=1`` (or ``--download-openmpi=1``).

.. note::

   Configuration of PETSc for a batch system requires special procedures described at the
   PETSc documentation site. One starts with a configure option ``--with-batch=1``. See
   the "Installing on machine requiring cross compiler or a job scheduler" section of the
   `PETSc installation page <PETSc-installation_>`_.

.. rubric:: Footnotes

.. [#f1] The ``PETSC_ARCH`` variable is just a string you can use to choose different
       PETSc configurations and does not have any other significance.
