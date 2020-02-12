.. include:: ../global.txt

.. _sec-install-spack:

Installing PISM using Spack
---------------------------

On supercomputers, Linux, and macOS PISM can be installed using the Spack_ package
manager.

Installing PISM using this method is easy: install Spack itself (see `Spack
documentation`_) and then run

.. code-block:: bash

   spack install pism

This will install PISM and all its prerequisites, including PETSc. The default PETSc
configuration in its Spack package includes many optional features not used by PISM. You
may want to disable these; to do this, use this command instead:

.. code-block:: bash

   spack install pism ^petsc~metis~hdf5~hypre~superlu-dist

.. note::

   The Spack package for PISM is maintained by Elizabeth Fischer (|efischer-email|);
   please e-mail her with questions about it. General questions about Spack_ should be
   sent directly to Spack developers.

.. note::

   With default Spack settings this installation method relies on building most (if not
   all) of PISM's prerequisites from sources.

   **This may take a very long time.**

   Please see `Spack documentation`_ (system packages) for a way to avoid re-building
   tools and libraries available on your system.
