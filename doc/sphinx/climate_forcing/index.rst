.. include:: shortcuts.txt

.. _sec-climate-forcing:

Climate forcing
===============

PISM has a well-defined separation of climate forcing from ice dynamics. This manual is
about the climate forcing interface.

By contrast, most options documented in the :ref:`sec-users-manual` control the ice
dynamics part. The User's Manual does, however, give an :ref:`overview of PISM's surface
(atmosphere) and ocean (sub-shelf) interfaces <sec-climate-inputs>`. At these interfaces
the surface mass and energy balances are determined and/or passed to the ice dynamics
code.

To get started with climate forcing usage we need to introduce some language to describe
parts of PISM. In this manual a *component* is a piece of PISM code, usually a C++ class.
A combination of components (or, in some cases, one component) makes up a "model" --- an
implementation of a physical/mathematical description of a system.

PISM's climate forcing code has two kinds of components.

- Ones that can be used as "stand-alone" models, such as the implementation of the PDD
  scheme (section :ref:`sec-surface-pdd`). These are *model components*.
- Ones implementing "corrections" of various kinds, such as lapse rate corrections
  (sections :ref:`sec-surface-lapse-rate` and :ref:`sec-atmosphere-lapse-rate`) or
  ice-core derived offsets (sections :ref:`sec-surface-delta-t` and
  :ref:`sec-ocean-delta-sl`, for example). These are called *modifier components* or
  *modifiers*.

Model components and modifiers can be chained as shown in
:numref:`fig-climate-input-data-flow`. For example,

.. code-block:: none

    -ocean constant,delta_T -ocean_delta_T_file delta_T.nc

combines the component providing constant (both in space and time) ocean boundary
conditions with a modifier that applies scalar temperature offsets. This combination
one of the many ocean models that can be chosen using components as building blocks.

Section :ref:`sec-forcing-examples` gives examples of combining components to choose
models. Before that we address how PISM handles model time (Section
:ref:`sec-model-time`).

.. admonition:: Summary of the main idea in using this manual

   Setting up PISM's climate interface *requires* selecting one surface and one ocean
   component. The surface component may use an atmosphere component also; see
   :numref:`fig-climate-input-data-flow`. Command-line options :opt:`-atmosphere`,
   :opt:`-surface` and :opt:`-ocean` each take a comma-separated list of keywords as an
   argument; the first keyword *has* to correspond to a model component, the rest can be
   "modifier" components. Any of these options can be omitted to use the default
   atmosphere, surface or ocean model components, but one has to explicitly choose a model
   component to use a modifier. Model components and modifiers are chained as in
   :numref:`fig-climate-input-data-flow`.

.. toctree::

   time.rst

   examples.rst

   testing.rst

   surface.rst

   atmosphere.rst

   ocean.rst
