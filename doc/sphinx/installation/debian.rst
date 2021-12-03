.. _sec-install-debian:

Installing prerequisites on Debian or Ubuntu
--------------------------------------------

You should be able to use your package manager to get the prerequisites for PISM. Install
the following packages using ``apt-get`` or ``synaptic`` or similar. All of these are
recommended as they satisfy requirements for building or running PISM.

.. csv-table:: Debian packages
   :name: tab-debian-packages
   :header: Name, Comment
   :widths: 2,5
   :file: debian-packages.csv

You may be able to install these by running

.. literalinclude:: code/install_libraries.sh
   :language: bash

.. only:: html

   Click :download:`here <code/install_libraries.sh>` to download this file.

(You may need to add ``sudo`` or change this command to match your package system.)

The command above takes care of all PISM prerequisites, including PETSc. Set
``PETSC_DIR=/usr/lib/petsc``\ [#petsc-arch]_ and follow the steps in
:ref:`sec-install-pism` to build PISM itself.

.. rubric:: Footnotes

.. [#petsc-arch] In this case you do not need to set ``PETSC_ARCH``.
