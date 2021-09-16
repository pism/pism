.. include:: ../../../global.txt

.. _sec-damage:

Fracture density (Continuum Damage Mechanics) modeling
------------------------------------------------------

.. contents::


The fracture-density approach in PISM based on :cite:`AlbrechtLevermann2012` and assumes
a macroscopic measure for the abundance of (partly microscale) crevasses and rifts that
form in ice (shelves) and that can be transported with the ice flow as represented in a
continuum ice-flow model. The approach is similar to the Continuum Damage Mechanics (CDM)
(e.g. :cite:`Lemaitre1996Damage` and :cite:`Borstad2013Creep`), where a damage state
variable `\phi` (or `D`) is introduced that equals zero for fully intact ice, and one for
fully fractured ice. This can be interpreted as a loss of all load bearing capacity.

The feedback of damage to the ice flow (creep) works within the existing constitutive
framework by introducing a linear mapping between the actual physical (damaged) state
of the material and an effective state that is compatible with a homogeneous,
continuum representation of the creep law (Eq. 6 in :cite:`AlbrechtLevermann2014softening`)

Fractures form above a critical stress threshold `\sigma_t` (or `\tau_0`) in the ice
(von Mises or maximum stress criterion or fracture toughness from Linear Elastic
Fracture Mechanics) with a certain fracture growth rate `\gamma`, that is related to the
strain rate (longitudinal spreading or effective strain rate;
Eq. 9 in :cite:`AlbrechtLevermann2012`). Fracture healing is assumed to occur with
a defined healing rate below a strain rate threshold (scaled with the difference to the
threshold or constant; Eq. 11 in :cite:`AlbrechtLevermann2012`). The fracture density
model in PISM thus requires four parameters (see option below). The fracture growth rate
`\gamma` is discarded if following the constitutive law by :cite:`Borstad2016Constitutive`.



.. list-table:: Options which turn on fracture model and softening effect.
   :name: tab-damage
   :header-rows: 1
   :widths: 1,3

   * - Option
     - Description

   * - :opt:`-fractures`
     - switch on fracture density model (fracture formation, transport and healing)
       according to :cite:`AlbrechtLevermann2012`

   * - :opt:`-fracture_parameters`
     - sets four parameters for fracture growth rate [], fracture stress threshold [Pa],
       fracture healing rate [] and strain rate healing threshold [s-1]

   * - :opt:`-constitutive_stress_limit`
     - uses fracture growth according to constitutive law in
       :cite:`Borstad2016Constitutive` (Eq. 4), discarding the choice of the
       fracture growth rate parameter

   * - :opt:`-constant_healing`
     - applies a constant healing rate `-\gamma_h \dot{\epsilon}_h` independent
       of local strain rate

   * - :opt:`-fracture_weighted_healing`
     - adds a term `1-D` to the healing term (similar to the source term),
       assuming that highly damaged ice heals slower, can be combined with :opt:`-constant_healing

   * - :opt:`-max_shear`
     - uses the maximum shear stress criterion for fracture formation (a.k.a. Tresca
       or Guest criterion in literature), which is more stringent than the default
       von-Mises criterion, see Eq. 7 in :cite:`AlbrechtLevermann2014softening`

   * - :opt:`-lefm`
     - uses the mixed-mode fracture toughness stress criterion based on Linear Elastic
       Fracture Mechanics, see Eqs. 8-9 in :cite:`AlbrechtLevermann2014softening`

   * - :opt:`-do_frac_on_grounded`
     - Model fracture density also in grounded ice areas (e.g. along ice stream shear zones)

   * - :opt:`-phi0`
     - Assuming that all ice entering the ice shelf from :var:`bc_mask` (e.g. via ice streams)
       has a constant fracture density

   * - :opt:`-scheme_fd2d`
     - Making use of advanced two-dimensional transport scheme to reduce the effect of
       numerical diffusion (Eq. 10 in :cite:`AlbrechtLevermann2014softening`)

  * - :opt:`-fracture_softening`
     - Parameter for the feedback strength of damage on the ice flow, if one no feedback,
       if zero full feedback (see `\epsilon` in Eq. 6 in :cite:`AlbrechtLevermann2014softening`)

  * - :opt:`-constant_fd`
     - Do not update the fracture density fields, but for instance make use of its softening effect.


In order to test for different damage options and parameters see the `example/ross/fracture` folder.
Build a setup for the Ross Ice Shelf and let the damage field evolve, with fracture bands reaching all 
the way from the inlets to the calving front.