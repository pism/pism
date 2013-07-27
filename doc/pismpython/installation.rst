.. _installation:

=============================================
Installation
=============================================

The Python bindings are an optional part of PISM, but 
are required to run the inversion code. The instructions
below show how to build a version of PISM with
support for Python bindings and inverse code assuming
that a working version of PISM has already been setup
(as described in the PISM Installation Manual
:cite:`pism-installation-manual`).

Prerequisites
=============

Python is already required to build PISM.  The following 
prerequisites will also be needed for the bindings
and inverse library.

``SWIG``
-----------

``SWIG`` (http://www.swig.org) is a software tool for automatically generating wrappers for C/C++ code. It can be be installed via

.. code-block:: bash

  sudo apt-get install swig

on Linux or 

.. code-block:: bash

  sudo port install swig

on Mac OS X.


``petsc4py``
------------

``petsc4py``  
is a set of Python bindings for PETSc.  `Download
<http://code.google.com/p/petsc4py/>`__
a version appropriate for your version of PETSc; version
3.3 works with PETSc 3.2 and 3.3.  Then, assuming you 
built PISM with ``PETSC_DIR=/home/user/petsc-3.2-p7/``
and ``PETSC_ARCH=linux-gnu-opt``

.. code-block:: bash

  export PETSC_DIR=/home/user/petsc-3.2-p7/
  export PETSC_ARCH=linux-gnu-opt
  python setup.py build
  python setup.py install --prefix INSTALL_LOCATION

This will install a copy of ``petsc4py`` in
:file:`{INSTALL_LOCATION}/lib/python{X}.{X}/site-packages`.
For example, for Python 2.7 and 
``INSTALL_LOCATION=/home/user/lib`` the install
location will be :file:`/home/user/lib/python2.7/site-packages`.
You will then want to add the install location to your ``PYTHONPATH``
by adding, e.g.,

.. code-block:: bash

  export PYTHONPATH=`/home/user/lib/python2.7/site-packages`:${PYTHONPATH}

to your :file:`.profile` or :file:`.bashrc`.

You can verify that you have installed ``petsc4py`` if the command

.. code-block:: bash

  python -c "import petsc4py"

Returns error free.

TAO
---

TAO is a software library, based on PETSc, for large-scale 
optimization problems.  Its use is optional; if it is not installed,
some of the inverse problem algorithms in PISM will simply not be available.

`Download <http://www.mcs.anl.gov/research/projects/tao/download/index.html>`_-
a copy of a version that is compatible with your version of PETSc 
(TAO 2.1 is compatible with PETSc 3.2 and 3.3). Then

.. code-block:: bash

  export PETSC_DIR=/home/user/petsc-3.2-p7/
  export PETSC_ARCH=linux-gnu-opt
  export TAO_DIR=`pwd`
  make all

You will need to add the ``TAO_DIR`` environment variable 
in your :file:`.profile` or :file:`.bashrc`, e.g.

.. code-block:: bash

  export TAO_DIR=/home/user/tao-2.1-p1

``Sphinx``
----------

Sphinx is a documentation generation tool and
can be installed via ``apt-get`` or ``macports``.
See the `installation instructions <http://sphinx-doc.org/latest/install.html>`_
for more details.  It is only required if you wish to 
build the python/inverse documentation.  For example, do

.. code-block:: bash

  sudo apt-get install sphinx-common

The documentation also requires the Sphinx extension called
``sphinxcontrib.bibtex``, which may come with some Sphinx packages (but not
with Debian packages at this time).  Without it you will see this error when
you try to build the documentation (see below):

.. code-block:: bash

  Extension error:
  Could not import extension sphinxcontrib.bibtex (exception: No module named bibtex)

To install it see the
`online instructions <http://sphinxcontrib-bibtex.readthedocs.org>`_.

Note that if you install Sphinx using macports,
you will install a version that depends on your python
version, and its executables will have names that
depend on the python version, e.g. ``sphinx-build-2.7``
rather than ``sphinx-build`` for Python 2.7.  You will want to
set up aliases so that the standard names work as well. To do this,

.. code-block:: bash

  sudo port select sphinx py27-sphinx

(replacing py27-sphinx with py26-sphinx for Python 2.6, etc.)

If you opt not to do this, you can tell ``cmake`` the
name of your sphinx executable using

.. code-block:: bash

  cmake -DSPHINX_EXECUTABLE=sphinx-build-2.7 ...


Building PISM with Python bindings
==================================

To setup a PISM build with Python bindings, either use

.. code-block:: bash

  cmake -DPism_BUILD_PYTHON_BINDINGS=1 ...

or, if using ``ccmake``, set ``Pism_BUILD_PYTHON_BINDINGS`` to ``ON``
in the user interface.

If ``cmake`` is unable to find ``petsc4py``, it will terminate
with the error 

.. code-block:: bash

  Could NOT find PETSc4Py (missing: PETSC4PY_INCLUDES)

If this occurs, verify that ``petsc4py`` can be found
in in your ``PYTHONPATH`` (i.e. ``python -c "import petsc4py"`` returns
error free).

PISM will link to TAO automatically if it can find a TAO installation.
If you wish to include TAO support and see the ``cmake`` log messages:

.. code-block:: bash

  -- Checking for package 'TAO'
  -- TAO_DIR is TAO_DIR-NOTFOUND

verify that you have set your ``TAO_DIR`` environment variable correctly.


Building the Documentation
==========================

In the PISM build directory, 

.. code-block:: bash

  make pismpython_docs

If you get an error like

.. code-block:: bash

  make: *** No rule to make target `pismpython_docs'.  Stop.

then re-run ``cmake ..`` or ``ccmake ..``, making sure that Sphinx is installed
(see above); the ``pismpython_docs`` target will then be present.
Once built, the main page for the documentation is then in
:file:`doc/pismpython/html/index.html` inside your build directory. The
documentation build can take some time while it
builds a large number of small images from
:math:`\text{\LaTeX}` formulas.

