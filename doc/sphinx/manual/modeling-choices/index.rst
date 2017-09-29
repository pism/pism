.. include:: ../../global.rst

.. _sec-modeling-choices:

Modeling choices
================

.. toctree::
   :maxdepth: 1

   computational/index.rst
   dynamics/index.rst
   subglacier/index.rst
   marine/index.rst



   
PISM is controlled by command-line options and configuration parameters. Every parameter
corresponds to an option; some options map to a combination of parameters.

There are five kinds of options and parameters:

#. Strings (mostly file names),
#. Scalars (a number representing a physical quantify, usually with units),
#. Integers (an integer number, usually corresponding to a count),
#. Booleans (True/False, on/off, yes/no),
#. Keywords (one of a given set of strings).

Each parameter can be set using the command-line option consisting of a dash followed by
the parameter name. For example,

.. code-block:: none

   -constants.standard_gravity 10

sets the acceleration due to gravity (parameter :config:`constants.standard_gravity`) to
`10`.

The following are equivalent and choose the temperature-based (as opposed to
enthalpy-based) energy balance model:

.. code-block:: none

   -energy.temperature_based
   -energy.temperature_based on
   -energy.temperature_based yes
   -energy.temperature_based true
   -energy.temperature_based True

The following are also equivalent: they disable updating geometry by performing a step of
the mass-continuity equation:

.. code-block:: none

   -geometry.update.enabled off
   -geometry.update.enabled no
   -geometry.update.enabled false
   -geometry.update.enabled False
   -no_geometry.update.enabled

The ``-no_`` prefix is still supported for compatibility with older scripts, but will be
removed in a later PISM version.
   
See FIXME for the full list of configuration parameters and corresponding options.
