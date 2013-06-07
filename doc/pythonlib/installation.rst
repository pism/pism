.. _installation:

=============================================
Python Bindings and Inverse Code Installation
=============================================

PISM's Python bindings can be used to write Python
scripts that use 

The Python bindings are an optional part of PISM, but 
are required to run the inversion code. The instructions
below show how to build a version of PISM with
support for Python bindings and inverse code assuming
that a working version of PISM has already been setup
(as described in the PISM Installation Manual
:cite:`pism-installation-manual`).

Python is already required to build PISM.  The following 
prerequisites will also be needed for the bindings
and inverse library.

``SWIG``
=========

``SWIG`` (http://www.swig.org) is a software tool for automatically generating wrappers for C/C++ code. It can be be installed via

.. code-block:: bash
  sudo apt-get install swig

on Linux or 

.. code-block:: bash
  sudo port install swig

on Mac OS X.

``Sphinx``
==========

Sphinx is a documentation generation tool and
can be installed via ``apt-get`` or ``macports``.
See the :ref:`installation instructions <http://sphinx-doc.org/latest/install.html>`
for more details.

``petsc4py``
============


TAO
===

