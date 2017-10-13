.. include:: ../global.txt

.. _sec-install-documentation:

Rebuilding PISM documentation
=============================

You might want to rebuild the documentation from source, as PISM and its
documentation evolve together. These tools are required:

.. list--table::
   :header-rows: 1

   * - Tool
     - Comment
   * - Sphinx_
     - needed for rebuilding this Manual
   * - ``sphinxcontrib.bibtex``
     - needed for rebuilding this Manual, as well as the documentation of PISM's Python
       bindings, below
   * - LaTeX_
     - needed for rebuilding the PDF version of this Manual and the source code browser
   * - ``latexmk``
     - needed for rebuilding the PDF version of this Manual
   * - doxygen_
     - required to rebuild the *Browser*
   * - graphviz_
     - required to rebuild the *Browser*

To install ``sphinxcontrib.bibtex`` see |sphinxcontrib-bibtex-url|.

Manual
------

Run the following commands to re-build the manual.

.. code-block:: bash

   # Change to the PISM source directory, then
   cd doc/sphinx
   make html     # to build the HTML manual
   make latexpdf # to build the PDF manual

The HTML manual will be in ``doc/sphinx/_build/html``. The PDF manual will be in
``doc/sphinx/_build/latex``.

Source Code Browser
-------------------

To rebuild the Source Code Browser, change to the PISM build directory and run

.. code-block:: bash

   make browser

To build documentation on a system without PISM’s prerequisite libraries (such as MPI and
PETSc), assuming that PISM sources are in ``~/pism-stable``, do the following:

.. code-block:: bash

   cd ~/pism-stable
   mkdir doc-build # create a build directory
   cd doc-build
   cmake ../doc

then commands "``make pism_manual``", "``make pism_forcing``" and others (see above) will
work as expected.

Documentation for PISM’s Python bindings and inversion tools
------------------------------------------------------------

The documentation for PISM’s Python bindings uses Sphinx_. The bindings make scripting and
interactive PISM possible, but many PISM users will not need them. Installing them is
required to use PISM for inversion of surface velocities for basal shear stress and ice
hardness. Building their documentation is strongly-recommended before use.

The bindings documentation also requires the Sphinx extension called
``sphinxcontrib.bibtex``, which may come with some Sphinx packages (but not with Debian
packages at this time). Without it you will see this error when you try to build the
bindings documentation:

.. code-block:: none

   Extension error:
   Could not import extension sphinxcontrib.bibtex (exception: No module named bibtex)


Note that if you install Sphinx using MacPorts_, you will install a version that depends
on your Python version, and its executables will have names that depend on the Python
version, e.g. ``sphinx-build-2.7`` rather than ``sphinx-build`` for Python 2.7. You will
want to set up aliases so that the standard names work as well. To do this,

.. code-block:: none

    sudo port select sphinx py27-sphinx

(replacing ``py27-sphinx`` with ``py26-sphinx`` for Python 2.6, etc.) If you opt not to do
this, you can tell CMake the name of your Sphinx executable using

.. code-block:: none

   cmake -DSPHINX_EXECUTABLE=sphinx-build-2.7 ...

for example.

Now you can build the documentation. In the PISM build directory, do

.. code-block:: none

    make pismpython_docs

If you get an error like

.. code-block:: none

   make: *** No rule to make target `pismpython_docs'.  Stop.

then re-run ``cmake ..`` or ``ccmake ..``, making sure that Sphinx is installed (see
above); the ``pismpython_docs`` make target will then be present.

The main page for the documentation is then in ``doc/pismpython/html/index.html`` inside
your build directory. The documentation build can take some time while it builds a large
number of small images from LaTeX formulas.
