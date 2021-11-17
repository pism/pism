.. include:: ../global.txt

.. _sec-install-pism:

Building PISM
-------------

To make sure that the key PETSc and MPI prerequisites work properly together, so that you
can run PISM in parallel, you might want to make sure that the correct ``mpiexec`` can be
found, by setting your ``PATH``. For instance, if you used the option
``--download-mpich=1`` in the PETSc configure, the MPI ``bin`` directory will have a path
like ``$PETSC_DIR/$PETSC_ARCH/bin``. Thus the following lines might appear in your
``.bashrc`` or ``.profile``, if not there already:

.. code-block:: bash

   export PETSC_DIR=/home/user/petsc-3.10.2/
   export PETSC_ARCH=opt
   export PATH=$PETSC_DIR/$PETSC_ARCH/bin/:$PATH

From now on we will assume that the ``PETSC_ARCH`` and ``PETSC_DIR`` variables are set.

.. note::

   The ``PETSC_ARCH`` variable is not needed if PETSc was configured using the
   ``--prefix=`` option.

Follow these steps to build PISM:

#. Get the latest source for PISM using the Git_ version control system by running

   .. _git-clone:

   .. code-block:: bash

      git clone git://github.com/pism/pism.git pism-stable

   A directory called "``pism-stable``" will be created. Note that in the future when you
   enter that directory, ``git pull`` will update to the latest revision of PISM. [#f1]_

   .. note::

      You can also `download a tarball from GitHub <pism-releases_>`_.

#. Build PISM:[#f2]_

   .. code-block:: bash

      mkdir -p pism-stable/build
      cd pism-stable/build
      export CC=mpicc
      export CXX=mpicxx
      cmake -DCMAKE_INSTALL_PREFIX=~/pism ..
      make -j install

   Here ``pism-stable`` is the directory containing PISM source code while ``~/pism`` is
   the directory PISM will be installed into.

   Variables ``CC`` and ``CXX`` specify MPI compiler wrappers provided by your MPI
   installation.

   .. note::

      When using MPI's compiler wrappers, make sure that ``mpicc`` and ``mpicxx`` you
      select were used to compile the PETSc library: *PISM and PETSc have to use the same
      MPI installation.*

   Commands above will configure PISM to be installed in ``~/pism/bin`` and
   ``~/pism/lib/`` then compile and install all its executables and scripts.

   If your operating system does not support shared libraries\ [#f3]_, then set
   ``Pism_LINK_STATICALLY`` to "ON". This can be done by either running

   .. code-block:: bash

      cmake -DPism_LINK_STATICALLY=ON ..

   or by using ``ccmake``\ [#f4]_ run

   .. code-block:: bash

      ccmake ..

   and then change ``Pism_LINK_STATICALLY`` (and then press ``c`` to "configure" and ``g``
   to "generate Makefiles"). Then run ``make install``.

   Temporary files created during the build process (located in the ``build``
   sub-directory) are not automatically deleted after installing PISM, so run "``make
   clean``" if space is an issue. You can also delete the build directory altogether if
   you are not planning on re-compiling PISM.

   .. note::

      When using Intel's compiler and high optimization settings such as ``-O3``,
      ``-fp-model precise`` may be needed to get reproducible model results. Set it using
      ``ccmake`` or by setting ``CFLAGS`` and ``CXXFLAGS`` environment variables when
      building PISM's prerequisites (such as PETSc) and PISM itself.

      .. code-block:: bash

         export CFLAGS="-fp-model precise"
         export CXXFLAGS="-fp-model precise"
         cmake [other options] ..

   .. note::

      To achieve best performance it can be useful to tell the compiler to target the
      "native" architecture. (This gives it permission to use CPU instructions that may
      not work on older CPUs.)

      .. code-block:: bash

         export CFLAGS="-march=native"
         export CXXFLAGS="-march=native"
         cmake [other options] ..

#. PISM executables can be run most easily by adding the ``bin/`` sub-directory in your
   selected install path (``~/pism/bin`` in the example above) to your ``PATH``. For
   instance, this command can be done in the Bash_ shell or in your ``.bashrc`` file:

   .. code-block:: bash

      export PATH=~/pism/bin:$PATH

#. Now see section :ref:`sec-install-quick-tests` or :ref:`sec-start` to continue.

.. _sec-install-pism-cmake-options:

PISM's build-time configuration
===============================

Some of PISM's features (the ones requiring additional libraries, for example) need to be
enabled when building PISM. This section lists important build-time options.

.. csv-table::
   :header: Option, Description

   ``CMAKE_BUILD_TYPE``, \"build type\": set to \"Debug\" for development
   ``BUILD_SHARED_LIBS``, build shared (as opposed to static) libraries (this is the default)
   ``Pism_LINK_STATICALLY``, set CMake flags to try to ensure that everything is linked statically
   ``Pism_LOOK_FOR_LIBRARIES``, specifies whether PISM should look for libraries (disable this on Crays)
   ``Pism_BUILD_EXTRA_EXECS``, build additional executables (needed to run ``make test``)
   ``Pism_BUILD_PYTHON_BINDINGS``, build PISM's Python bindingd; requires ``petsc4py``
   ``Pism_USE_PROJ``, use the PROJ_ library to compute latitudes and longitudes of grid points
   ``Pism_USE_PIO``, use the ParallelIO_ library to write output files
   ``Pism_USE_PARALLEL_NETCDF4``, use NetCDF_ for parallel file I/O
   ``Pism_USE_PNETCDF``, use PnetCDF_ for parallel file I/O
   ``Pism_DEBUG``, enables extra sanity checks in the code (this makes PISM a lot slower but simplifies development)

To enable PISM's use of PROJ_, for example, run

.. code-block:: bash

   cmake -DPism_USE_PROJ [other options] ..

.. _sec-install-local-libraries:

Building PISM with libraries in non-standard locations
======================================================

To build PISM with libraries installed in a non-standard location such as ``~/local/``,
use CMake's variable ``CMAKE_FIND_ROOT_PATH``. Set it to a semicolon-separated list of
directories.

For example, if ``netcdf.h`` is located in ``~/local/netcdf/include/`` and
``libnetcdf.so`` is in ``~/local/netcdf/lib``, add ``~/local/netcdf`` to
``CMAKE_FIND_ROOT_PATH``:

.. code-block:: bash

   cmake -DCMAKE_FIND_ROOT_PATH=~/local/netcdf [other options] ..

To build PISM using parallel I/O libraries installed as described in
:ref:`sec-install-parallel-io-libs`, do this:

.. code-block:: bash

   cmake -DCMAKE_FIND_ROOT_PATH="~/local/netcdf;~/local/pnetcdf;~/local/parallelio" \
         -DPism_USE_PNETCDF \
         -DPism_USE_PARALLEL_NETCDF4 \
         -DPism_USE_PIO \
         ..

.. rubric:: Footnotes

.. [#f1] Of course, after ``git pull`` you will ``make -C build install`` to recompile and
       re-install PISM.

.. [#f2] Please report any problems you meet at these build stages by `sending us
       <pism-email_>`_ the output.

.. [#f3] This might be necessary if you’re building on a Cray XT5 or a Sun Opteron Cluster,
       for example.

.. [#f4] Install the ``cmake-curses-gui`` package to get ``ccmake`` on Ubuntu_.
