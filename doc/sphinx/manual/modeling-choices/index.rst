.. include:: ../../global.txt

.. _sec-modeling-choices:

Modeling choices
================

PISM consists of several sub-models corresponding to various physical processes and
parameterizations. In short, these are

#. Ice dynamics and thermodynamics (stress balance, ice rheology, mass and energy
   conservation, ice age)
#. Subglacial processes (hydrology, basal strength, bed deformation)
#. Marine ice-sheet modeling (parameterization of calving processes, calving front advance
   and retreat, iceberg removal)

All these sub-models are controlled by command-line options and configuration parameters.\
[#]_

In additions to this, one has to choose the computational grid, the modeling domain, and
so on.

This section describes how to understand and make these modeling choices.

.. toctree::
   :maxdepth: 1

   computational/index.rst

   dynamics/index.rst

   subglacier/index.rst

   marine/index.rst

   regional/index.rst

   disabling.rst

   hard-choices.rst

.. [#] See :ref:`sec-parameter-list` for the full list of configuration parameters and
       corresponding options.

