.. include:: ../global.txt

.. default-role:: literal

.. _sec-development-workflow:

Development workflow
====================

The recommended development workflow is:

#. When starting a new project, create a topic branch starting from `dev`. If you are
   fixing a bug in a released version of PISM, create a topic branch starting from
   `master`.
#. Make changes to the code or documentation (or both).

   a) Compile.
   b) Fix any compilation errors and warnings. Repeat until your code compiles without
      warnings.
   c) Push your code to GitHub for review or to get help with b).

#. Add verification or regression tests.
#. Test your code and repeat 2a and 2b until all tests pass.

#. Update documentation (if necessary).
#. Update the change log `CHANGES.rst`.\ [#]_
#. Merge new features into `dev` and fixes into `master` and `dev` (or submit a pull
   request).

This document covers the tools and approaches we found useful for the steps listed above.

.. contents::

.. _sec-editor:

Editing source code
^^^^^^^^^^^^^^^^^^^

Any text editor supporting C++ and Python will work, but we recommend using Emacs_.

Your editor needs to provide the ability to jump from a compiler's error message to the
relevant part of the code. In Emacs, use `M-x compile` to start a compilation and `M-x
recompile` to re-run it.

.. _sec-developmental-environment:

Setting up the environment for PISM development
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Configure PISM with debugging symbols enabled.

.. code-block:: bash

   export PETSC_DIR=~/local/petsc/petsc-3.11.3/
   export PETSC_ARCH=opt

   CC=mpicc CXX=mpicxx cmake \
       -DCMAKE_BUILD_TYPE=Debug \
       -DPism_BUILD_EXTRA_EXECS=YES \
       -DPism_BUILD_PYTHON_BINDINGS=YES \
       -DPism_DEBUG=YES \
       ${pism_source_dir}

.. list-table:: PISM's configuration flags for development
   :name: tab-cmake-development-flags
   :header-rows: 1

   * - Flag
     - Meaning
   * - `-DCMAKE_BUILD_TYPE=Debug`
     - Enables pedantic compiler warnings
   * - `-DPism_BUILD_EXTRA_EXECS=YES`
     - Build extra testing executables (needed by some of regression test)
   * - `-DPism_BUILD_PYTHON_BINDINGS=YES`
     - Build PISM's Python bindings (used by many regression tests)
   * - `-DPism_DEBUG=YES`
     - Enables extra sanity checks in PISM

In most cases it is not necessary to build PETSc with debugging enabled.

.. _sec-compiling-pism:

Compiling PISM
^^^^^^^^^^^^^^

If the computer you use for development has multiple CPU cores you should tell `make` to
use all of them. Run `make -j4` on a four-core laptop, for example; this will
significantly speed up compilation.

To further speed up re-compiling PISM, install ccache_ and configure PISM as follows:

.. code-block:: bash

   CC="ccache mpicc" CXX="ccache mpicxx" cmake ...

.. _sec-debugging-pism:

Debugging
^^^^^^^^^

A debugger such as GDB_ or LLDB_ can be very useful.\ [#]_ There are many online tutorials
for both.

Some of the more troublesome bugs involve memory access errors (*segmentation fault*
errors are often caused by these). Consider using Valgrind_ to detect them.

.. note::

   Your code will run much, much slower when using Valgrind, so it is important to find
   a small test case reproducing the error.

.. _sec-continuous-integration:

Debugging issues caught by automatic tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every time somebody pushes changes to PISM's repository on GitHub the continuous
integration system attempts to build PISM and (if it was built successfully) run a suite
of tests.

It is often helpful to be able to run the same tests locally. To do this, install Docker_
and `CircleCI CLI`_ (command-line interface), then run

.. code-block:: none

   circleci local execute

in PISM's source code directory.

.. _sec-editing-the-manual:

Editing PISM's manual
^^^^^^^^^^^^^^^^^^^^^

PISM's manual is written using the reStructuredText_ markup and Sphinx_.

See :ref:`sec-install-documentation` for a list of tools needed to build PISM's
documentation.

When working on major edits, sphinx-autobuild_ can save you a lot of
time. Run

.. code-block:: bash

   sphinx-autobuild -B /path/to/pism/doc/sphinx/ /tmp/pism-manual

To get a browser window containing PISM's manual that will stay up to date with your
edits.

Run `make manual_html` to re-build the manual and open `doc/sphinx/html/index.html`
inside your build directory.

Lists of configuration parameters and diagnostics reported by PISM are generated from the
|config-cdl| and the PISM code itself. To make sure that they are up to date, run `make`
in `doc/sphinx` and commit the changes.

.. _sec-writing-tests:

Writing tests
^^^^^^^^^^^^^

All contributions containing new features should contain tests for the new code.\ [#]_

A contribution fixing a bug should (ideally) contain a test that will ensure that it is
fixed.

Add verification tests (tests comparing results to an analytical solution) whenever
possible. If a verification test is not an option, consider adding a *regression* test
that compares computed results to a stored output from a trusted version of the code. This
will make it easier to detect a regression, i.e. an undesirable change in model results.

Python bindings make it possible to test many PISM's components in isolation from the rest
of the code. See tests in `test/regression` for some examples.

.. note::

   This manual should cover PISM's Python bindings. If you see this, please e-mail
   |pism-email| and remind us to document them.

.. rubric:: Footnotes

.. [#] See `Keep a change log <keep-a-change-log_>`_ for inspiration.
.. [#] In most cases a serial debugger is sufficient.
.. [#] Contributions of tests for existing code are also welcome.
