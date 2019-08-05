.. _sec-install-debian:

Installing prerequisites from Debian packages
---------------------------------------------

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
   :lines: 3-

.. only:: html

   Click :download:`here <code/install_libraries.sh>` to download this file.

(You may need to change this command to match your package system.)

Once done, see :ref:`sec-install-petsc` to install PETSc from source and then
:ref:`sec-install-pism` for building PISM itself.
