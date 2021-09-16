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
   c) Run `make test` and fix any test failures.
   d) Push your code to GitHub for review or to get help with b) and c).

#. Add verification or regression tests.
#. Test your code and repeat 2a--2c until all tests pass.

#. Update documentation (if necessary).
#. Update the change log `CHANGES.rst`.\ [#]_
#. Merge new features into `dev` and fixes into `master` and `dev` (or submit a pull
   request).

This document covers the tools and approaches we found useful for the steps listed above.

.. contents::

.. _sec-developmental-environment:

Setting up the environment
^^^^^^^^^^^^^^^^^^^^^^^^^^

The majority of interesting PISM runs are performed on supercomputers, but we **do not**
recommend using supercomputers for development.

  Use a desktop (or a laptop) computer running Linux or macOS.

While you can use `SSH` to connect to a remote system to write, compile, and test your
code, doing so will reduce your productivity when compared to using a computer you have
physical access to.

Any MPI implementation would work, but we prefer to use MPICH_ for PISM development. This
MPI implementation

- has headers that compile without warnings,
- provides type checking for pointer arguments in MPI calls, and
- does not produce "false positives" when debugging memory access with Valgrind_.

When working on a fix for a stubborn bug it *may* be helpful to use PETSc compiled with
debugging enabled (option `--with-debugging=1`), but in our experience this is rarely
needed. Optimized PETSc builds (using `--with-debugging=0`) are faster and this helps with
overall productivity.

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

.. _sec-editor:

Editing source code
^^^^^^^^^^^^^^^^^^^

Any text editor supporting C++, Python, and reStructuredText will work, but we recommend
Emacs_.

Your editor needs to provide the ability to jump from a compiler's error message to the
relevant part of the code. In Emacs, use `M-x compile` to start a compilation and `M-x
recompile` to re-run it.

An editor that can help you navigate the code and find function definitions, etc is also
helpful; try an IDE such as KDevelop_, for example.

.. _sec-compiling-pism:

Compiling PISM
^^^^^^^^^^^^^^

If the computer you use for development has multiple CPU cores you should tell `make` to
use all of them. Run `make -j4` on a four-core laptop, for example; this will
significantly speed up compilation.

To further speed up re-compiling PISM, install ccache_ and configure PISM as follows:

.. code-block:: bash

   CC="ccache mpicc" CXX="ccache mpicxx" cmake ...

It may be helpful to use LLD_ to link PISM during development since it is a lot faster
than GNU ld. Add the following CMake_ options to give this a try.

.. code-block:: bash

   -DCMAKE_EXE_LINKER_FLAGS="-fuse-ld=lld" \
   -DCMAKE_SHARED_LINKER_FLAGS="-fuse-ld=lld" \
   -DCMAKE_MODULE_LINKER_FLAGS="-fuse-ld=lld"


.. _sec-debugging-pism:

Debugging
^^^^^^^^^

The first step in debugging an issue is *always* this:

  find the shortest simulation using the smallest possible grid that exhibits the
  problematic behavior.

It does not have to be *the* shortest simulation, but it should complete (or stop because
of a failure) within seconds when running on the machine used for development.

A debugger such as GDB_ or LLDB_ can be very useful.\ [#]_ There are many online tutorials
for both.

You will need to know how to

- start a program,
- interrupt execution,
- set and remove a breakpoint,
- continue execution after stopping at a breakpoint,
- continue execution to the next line of the code,
- continue execution to the end of the current function call,
- step into a function call,
- print the value of a variable,
- print the stack trace.

This basic set of debugging skills is often sufficient.

Sometimes a failure happens in a loop that iterates over grid points and stepping through
the code in a debugger is impractical. A *conditional breakpoint* would help (i.e. stop
only if a condition is true), but this debugger feature is not always well supported and
often significantly slows down execution.

Here's a different way to stop the code when a condition is met: add `#include <cassert>`
to the top of the file (if it is not there), then add `assert(!condition);` to the place
in the code where you would like to stop if `condition` is met.

For example,

.. code-block:: c++

   assert(!(i == 228 and j == 146));

will stop execution at the grid point where `i == 228` and `j == 146`.

Some of the more troublesome bugs involve memory access errors (*segmentation fault*
errors are often caused by these). Consider using Valgrind_ to detect them.

.. note::

   Your code will run much, much slower when using Valgrind, so it is important to find
   a small test case reproducing the error.

Issues visible in parallel runs only
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every once in a while a bug shows up in a parallel run but *not* in an equivalent serial
one. These bugs tend to be hard to fix and there is no definitive technique (or tool) that
helps with this. Here are some tips, though.

- Reduce the number of processes as much as possible. Most of the time the number of
  processes can be reduced all the way down to 2 (the smallest truly parallel case).
- Run PISM with the option `-start_in_debugger`. This will produce a number of terminal
  windows with GDB_. You will need to *continue execution* (GDB's command `c`) in all of
  the windows. If PISM freezes, interrupting execution and printing the stack trace would
  tell you where it got stuck.

  Executing commands in all the windows with GDB is tedious and error-prone. To execute a
  number of commands in all of them at the beginning of the run, create a file called
  `.gdbinit` (in the current directory) and put GDB commands there (one per line).

  For example,

  .. code-block:: none

     break pism::RuntimeError::RuntimeError()
     continue

  will set a breakpoint at `pism::RuntimeError::RuntimeError()` and continue execution.
- A parallel debugger such as TotalView may be helpful but requires a license. We don't
  have experience with it and cannot give any advice.

.. _sec-continuous-integration:

Issues caught by automatic tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every time somebody pushes changes to PISM's repository on GitHub the continuous
integration system attempts to build PISM and (if it was built successfully) run a suite
of tests.

It is often helpful to be able to run the same tests locally. To do this, install Docker_
and `CircleCI CLI`_ (command-line interface), then run

.. code-block:: bash

   circleci local execute --job={job}
   # where job is one of
   # build-gcc build-clang
   # build-clang-minimal build-gcc-minimal build-manual

in PISM's source code directory.

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

Here are some test writing tips:

- Make sure that a verification test uses a grid that is not square (different number of
  grid points in `x` and `y` directions).
- If possible, write a test that performs a number of computations along a *refinement
  path* and compare the computed convergence rate to the theoretical one.
- Try to include tests ensuring that `x` and `y` directions are interchangeable: in most
  cases flow from left to right should behave the save as from bottom towards the top, etc.

  Here are two ways to do this:

  - Repeat the test twice, the second time using transposed inputs, then transpose results
    and compare.
  - Repeat the test twice, once using a refinement path along `x` and the second time
    along `y`; make sure that you see the same convergence rate.
- It is good to check if the implementation preserves symmetries if the setup has any.
- If a test uses a temporary file, make sure that it will not clash with names of files
  used by other tests. One easy way to do this is by generating a unique file name using
  ``mktemp`` (in Bash scripts) or ``str(uuid.uuid4())`` (in Python).

Python bindings make it possible to test many PISM's components in isolation from the rest
of the code. See tests in `test/regression` for some examples.

.. note::

   This manual should cover PISM's Python bindings. If you see this, please e-mail
   |pism-email| and remind us to document them.

Running tests
~~~~~~~~~~~~~

Run ``make test`` in parallel by adding

.. code-block:: bash

   export CTEST_PARALLEL_LEVEL=N

to your ``.bashrc``. This will tell ``ctest`` to run ``N`` at the same time. Or run
``ctest -j N`` instead of ``make test``.

.. _sec-editing-the-manual:

Editing PISM's manual
^^^^^^^^^^^^^^^^^^^^^

PISM's manual is written using the reStructuredText_ markup and Sphinx_.

See :ref:`sec-install-documentation` for a list of tools needed to build PISM's
documentation.

When working on major edits, sphinx-autobuild_ can save you a lot of
time. Run

.. code-block:: bash

   make manual_autobuild

in the build directory to get a browser window containing PISM's manual that will stay up
to date with your edits.

To generate the HTML version of the manual, run `make manual_html`; the output is saved to
`doc/sphinx/html/` in your build directory.

Edit ``doc/sphinx/math-definitions.tex`` to add custom LaTeX commands used in formulas
(one per line).

.. _sec-manual-pism-parameters:

Listing configuration parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The list in :ref:`sec-parameter-list` is generated from |config-cdl|. Edit this file to
update the list.

When documenting a sub-model, use the ``:config:`` role to mention a parameter. This will
create a hyperlink to the complete list of parameters and ensure that all parameter names
are spelled correctly.

To create a list of parameters controlling a sub-model, use the ``pism-parameters``
directive. For example, to list all parameters with the prefix ``constants.ice.``, add
this:

.. code-block:: rst

   .. pism-parameters::
      :prefix: constants.ice.

This is the resulting list:

.. pism-parameters::
   :prefix: constants.ice.

For this to work, configuration parameters should be documented in |config-cdl| (see the
``..._doc`` attribute for each configuration parameter).

.. _sec-manual-pism-diagnostics:

Listing diagnostic quantities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The list of diagnostics reported by PISM is generated by running the PISM code itself. To
make sure that it is up to date, run `make` in `doc/sphinx` and commit the changes.

.. rubric:: Footnotes

.. [#] See `Keep a change log <keep-a-change-log_>`_ for inspiration.
.. [#] In most cases a serial debugger is sufficient.
.. [#] Contributions of tests for existing code are also welcome.
