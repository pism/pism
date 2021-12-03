.. include:: ../global.txt

.. _sec-install-documentation:

Rebuilding PISM documentation
=============================

You might want to rebuild the documentation from source, as PISM and its documentation
evolve together. These tools are required:

.. list-table::
   :header-rows: 1

   * - Tool
     - Comment
   * - Sphinx_ version 3.0 or newer
     - needed to rebuild this Manual
   * - sphinxcontrib.bibtex_
     - needed to rebuild this Manual
   * - LaTeX_
     - needed to rebuild the PDF version of this Manual
   * - Latexmk_
     - needed to rebuild the PDF version of this Manual
   * - ``ncgen`` from NetCDF_
     - needed to rebuild this Manual
   * - netcdf4-python_
     - needed to rebuild this Manual
   * - doxygen_
     - required to rebuild the |pism-browser|
   * - graphviz_
     - required to rebuild the |pism-browser|

On a Debian-based system you may be able to install most of these by running

.. literalinclude:: code/install_docu_libraries.sh
   :language: bash
   :lines: 3-

.. only:: html

   Click :download:`here <code/install_docu_libraries.sh>` to download this file.

(You may need to change this command to match your package system.)

Manual
------

To rebuild this manual, change to the PISM build directory and run

.. code-block:: bash

   make manual_html # for the HTML version of the manual
   make manual_pdf  # for the PDF version of the manual

The main page for this manual is then in ``doc/sphinx/html/index.html`` inside your build
directory.

The PDF manual will be in ``doc/sphinx/pism_manual.pdf`` in your build directory.

Source Code Browser
-------------------

To rebuild the |pism-browser|, change to the PISM build directory and run

.. code-block:: bash

   make browser

The main page for the documentation is then in ``doc/browser/html/index.html`` inside
your build directory.

Re-building documentation without PISM's prerequisites
------------------------------------------------------

To build documentation on a system without PISM’s prerequisite libraries (such as MPI and
PETSc), assuming that PISM sources are in ``~/pism-stable``, do the following:

.. code-block:: bash

   cd ~/pism-stable
   mkdir doc-build # create a build directory
   cd doc-build
   cmake ../doc

then commands "``make manual_html``", "``make manual_pdf``" and others (see above) will
work as expected.
