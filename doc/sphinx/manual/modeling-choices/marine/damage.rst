.. include:: ../../../global.txt

.. _sec-damage:

Fracture density modeling
-------------------------

The fracture density approach in PISM is based on :cite:`AlbrechtLevermann2012`
and assumes a macroscopic measure for the abundance of (partly microscale) crevasses
and rifts that form in ice (shelves) and that can be transported with the ice flow
as represented in a continuum ice-flow model. This approach is similar to the
Continuum Damage Mechanics (CDM) (e.g. :cite:`Lemaitre1996Damage` and
:cite:`Borstad2013Creep`) introducing a damage state variable (`\phi` or `D`) that
equals zero for fully intact ice and one for fully fractured ice, that can be
interpreted as a loss of all load bearing capacity.

The feedback of damage to the ice flow (creep) works within the existing constitutive
framework by introducing a linear mapping between the actual physical (damaged) state
of the material and an effective state that is compatible with a homogeneous,
continuum representation of the creep law (Eq. 6 in :cite:`AlbrechtLevermann2014softening`).

Fractures form above a critical stress threshold `\sigma_{\text{cr}}` in the ice (e.g. von
Mises criterion, maximum stress criterion or fracture toughness from Linear Elastic Fracture
Mechanics) with a fracture growth rate proportional to `\gamma` (Eq. 2 in
:cite:`AlbrechtLevermann2014softening`), that is related to the strain rate (longitudinal
spreading or effective strain rate; Eq. 9 in :cite:`AlbrechtLevermann2012`). Fracture
healing is assumed to occur with a defined healing rate below a strain rate threshold
(scaled with the difference to the threshold or constant; Eq. 11 in
:cite:`AlbrechtLevermann2012`).

The fracture growth constant `\gamma` (:config:`fracture_density.gamma`) is ignored if
:config:`fracture_density.borstad_limit` is set.

To enable this model, set :config:`fracture_density.enabled`.

.. rubric:: Parameters

Prefix: ``fracture_density.``

.. pism-parameters::
   :prefix: fracture_density.
   :exclude: fracture_density.enabled

.. rubric:: Testing

See the scripts in ``example/ross/fracture`` for a way to test different damage options
and parameter values. Build a setup for the Ross Ice Shelf and let the damage field
evolve, with fracture bands reaching all the way from the inlets to the calving front.
