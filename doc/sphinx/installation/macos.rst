.. include:: ../global.txt

.. _sec-install-macos:

Installing required libraries on macOS
--------------------------------------

Follow these steps to install PISM's prerequisites on the macOS operating system.

#. As PISM is distributed as source code only, you will need software developer’s tools,
   XCode_ and the *X window system server*, XQuartz_.

#. The use of MacPorts_ (or Fink_, or Homebrew_) is recommended, as it significantly
   simplifies installing many open-source libraries. These instructions assume that you
   use MacPorts_. Download a package from the MacPorts_, install, and set the environment:

   .. code-block:: bash

      export PATH=/opt/local/bin:/opt/local/sbin:$PATH

#. It may not be necessary to install Python, as it is bundled with the operating system.
   Some PISM scripts use SciPy; it can be installed using MacPorts or by downloading the
   `Enthought Python Distribution <Enthought_>`_.

#. This MacPorts command should install all of PISM's required libraries:

   .. code-block:: bash

      sudo port install git cmake fftw-3 gsl mpich-default netcdf udunits2 libproj4 ncview

#. At this point, all the PISM prerequisites except PETSc are installed. Proceed to
   :ref:`sec-install-petsc`.

#. Now you can build PISM as described in section :ref:`sec-install-pism`.
