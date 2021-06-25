.. include:: ../global.txt

.. _sec-install-quick-tests:

Quick tests of the installation
===============================

Once you’re done with the installation, a few tests can confirm that PISM is functioning
correctly.

#. Try a MPI four process verification run:

   .. code-block:: bash

      mpiexec -n 4 pismv -test G -y 200

   If you see some output and a final ``Writing model state`` ``to file ’unnamed.nc’``
   then PISM completed successfully. At the end of this run you get measurements of the
   difference between the numerical result and the exact solution. See :ref:`sec-verif`
   for more on PISM verification.

   The above "``-n 4``" run should work even if there is only one actual processor (core)
   on your machine. (In that case MPI will just run multiple processes on the one
   processor.) This run will also produce a NetCDF output file ``unnamed.nc``, which can
   be read and viewed by NetCDF tools.

#. Try an EISMINT II run using the PETSc viewers (under the X window system):

   .. code-block:: none

      pismr -eisII A -y 5000 -view thk,temppabase,velsurf_mag

   When using such viewers and ``mpiexec`` the additional final option ``-display :0`` is
   sometimes required to enable MPI to use X, like this:

   .. code-block:: none

       mpiexec -n 2 pismr -eisII A -y 5000 -view thk,temppabase,velsurf_mag -display :0

   Also ``-drawpause 0.1`` or similar may be needed if the figures are refreshing too fast.

#. Run a basic suite of software tests. To do this, make sure that NCO_ and Python
   packages NumPy_ and netcdf4-python_ are installed. Also, the CMake flag
   ``Pism_BUILD_EXTRA_EXECS`` should be ``ON``. Then run:

   .. code-block:: bash

      make       # do this if you changed something with CMake
      make test

   in the build directory.

   The message at the bottom of the output should say

      ``100% tests passed, 0 tests failed out of XX``

   or similar.

   Feel free to `e-mail us <pism-email_>`_ about any test failures you see. Please run

   .. code-block:: bash

      ctest --output-on-failure > make-test.log

   and send us the ``make-test.log`` that this produces.

Next steps
==========

Start with the section :ref:`sec-start`.

Completely up-to-date documentation can be built from source; see
:ref:`sec-install-documentation` for details.

A final reminder with respect to installation: Let’s assume you have checked out a copy of
PISM using Git, :ref:`as described above <git-clone>`. You can then update your copy of PISM
to the latest version by running ``git pull`` in the PISM directory and ``make install``
in your build directory.
