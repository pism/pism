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

Follow these steps to build PISM:

#. Get the latest source for PISM using the Git_ version control system:

   Check `PISM's website <PISM_>`_ for the latest version of PISM.

   .. _git-clone:

   Run

   .. code-block:: bash

      git clone git://github.com/pism/pism.git pism-stable

   A directory called "``pism-stable``" will be created. Note that in the future when you
   enter that directory, ``git pull`` will update to the latest revision of PISM. [#]_

#. Build PISM:[#]_

   .. code-block:: bash

      mkdir -p pism-stable/build
      cd pism-stable/build
      PISM_INSTALL_PREFIX=~/pism CC=mpicc CXX=mpicxx cmake ..
      make install

   Here ``pism-stable`` is the directory containing PISM source code while ``~/pism`` is
   the directory PISM will be installed into.

   Variables ``CC`` and ``CXX`` specify MPI compiler wrappers provided by your MPI
   installation.

   .. note::

      When using MPI's compiler wrappers, make sure that ``mpicc`` and ``mpicxx`` you
      select were used to compile the PETSc library: *PISM and PETSc have to use the same
      MPI installation.*

   All the temporary files created during the build process will be in
   ``pism-stable/build`` created above.

   Commands above will configure PISM to be installed in ``~/pism/bin`` and
   ``~/pism/lib/`` then compile and install all its executables and scripts.

   If your operating system does not support shared libraries\ [#]_, then set
   ``Pism_LINK_STATICALLY`` to "ON". This can be done by either running

   .. code-block:: bash

      cmake -DPism_LINK_STATICALLY=ON ..

   or by using ``ccmake``\ [#]_ run

   .. code-block:: bash

      ccmake ..

   and then change ``Pism_LINK_STATICALLY`` (and then press ``c`` to "configure" and ``g``
   to "generate Makefiles"). Then run ``make install``.

   Object files created during the build process (located in the ``build`` sub-directory)
   are not automatically deleted after installing PISM, so run "``make clean``" if space
   is an issue. You can also delete the build directory altogether if you are not planning
   on re-compiling PISM.

   .. note::

      When using Intel's compiler high optimization settings such as ``-O3``, ``-fp-model
      precise`` may be needed to get reproducible model results. Set it using ``ccmake``
      or by setting ``CFLAGS`` and ``CXXFLAGS`` environment variables when building PISM's
      prerequisites and PISM itself.

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


.. rubric:: Footnotes

.. [#] Of course, after ``git pull`` you will ``make -C build install`` to recompile and
       re-install PISM.

.. [#] Please report any problems you meet at these build stages by `sending us
       <pism-email_>`_ the output.

.. [#] This might be necessary if you’re building on a Cray XT5 or a Sun Opteron Cluster,
       for example.

.. [#] Install the ``cmake-curses-gui`` package to get ``ccmake`` on Ubuntu_.
