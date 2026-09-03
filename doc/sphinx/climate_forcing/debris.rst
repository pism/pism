.. include:: shortcuts.txt

.. _sec-debris:

Supraglacial debris
-------------------

.. contents::

PISM's debris models provide the thickness of the supraglacial debris layer
:var:`debris_thickness` and, through the melt enhancement component (see
:ref:`sec-debris-melt-enhancement`), its effect on ice melt. They are selected with
:opt:`-debris` (:config:`debris.models`), a comma-separated list of a model and any number
of modifiers, as for the other climate forcing components. Leave the list empty (the
default) to disable debris modeling.

.. note::

   In the current version the melt enhancement factor is *computed and reported* but not
   yet applied by the surface mass balance models, i.e. the debris is a passive tracer.

.. _sec-debris-given:

Reading the debris thickness from a file
++++++++++++++++++++++++++++++++++++++++

:|options|: ``-debris given``
:|variables|: :var:`debris_thickness` (m)
:|implementation|: ``pism::debris::Given``

Reads the (possibly time-dependent) debris thickness from :config:`debris.given.file`
(:opt:`-debris_given_file`); set :config:`debris.given.periodic` to treat the forcing as
periodic in time.

.. _sec-debris-transport:

Prognostic debris transport
+++++++++++++++++++++++++++

:|options|: ``-debris transport``
:|variables|: :var:`debris_input_rate` (m of solid debris per year, in
              :config:`debris.transport.input.file`)
:|implementation|: ``pism::debris::DebrisTransport``

This model follows :cite:`VerhaegenHuybrechts2026`. Debris supplied by the surrounding
terrain at the prescribed rate :var:`debris_input_rate` is

- buried by accumulation where the surface mass balance is positive and advected
  englacially by the three-dimensional ice velocity (equations 15 and 16 in the
  paper),
- released to the surface where ice melts (melt-out, equation 19), or added to the
  surface directly where the mass balance is negative (equation 18),
- advected by the ice surface velocity, redistributed down the slope of the debris
  surface `h_s + h_d` by the gravitational flux `F_G = -K_d \nabla(h_s + h_d)` with
  `K_d = \mu_d (1 - \phi_d) \rho_d g h_d` (equations 21 and 22), and
- removed at glacier-margin cells sloping into the ice-free foreland at the rate
  `|F_G| / \Gamma`, where `\Gamma` is the marginal length scale (equation 23). If
  `\Gamma` exceeds the grid spacing, the debris thickness and the slope are averaged over
  the margin cell and `\lceil \Gamma / \Delta x \rceil - 1` cells up-glacier along the
  steepest ascent (at most two).

The englacial debris is stored as mass per unit area in control volumes around the
levels of PISM's vertical grid and advected in flux form, so the total debris mass is
conserved to rounding error; ice removed at the surface, at the base, or in cells that
become ice-free releases or removes the debris it contains. Unlike the paper, PISM does
not rescale the debris mass to enforce conservation; instead it reports a mass budget (see
below). The transport schemes are selected with :config:`debris.transport.englacial.scheme`
(``upwind`` or ``mpdata``) and :config:`debris.transport.supraglacial.scheme` (``upwind``,
``mpdata``, ``uno2``, ``uno3``); MPDATA :cite:`Smolarkiewicz1983` is the anti-diffusive
scheme used in the paper, with :config:`debris.transport.mpdata.iterations` passes and the
non-oscillatory limiter enabled by :config:`debris.transport.mpdata.nonoscillatory`.

The explicit vertical advection restricts the time step to a fraction
(:config:`debris.transport.vertical_cfl_ratio`) of the smallest vertical spacing divided
by the largest vertical velocity; the gravitational transport is sub-cycled internally
(:config:`debris.transport.diffusion_cfl_ratio`).

.. list-table:: Parameters of the transport model
   :name: tab-debris-transport
   :header-rows: 1
   :widths: 1,1

   * - Parameter
     - Meaning
   * - :config:`debris.transport.density`, :config:`debris.transport.porosity`
     - Density `\rho_d` and porosity `\phi_d` of the debris.
   * - :config:`debris.transport.mobility`
     - Efficiency `\mu_d` of the gravitational transport.
   * - :config:`debris.transport.marginal_length_scale`
     - `\Gamma`, the length scale of the removal at the margin.
   * - :config:`debris.transport.gravitational_transport`,
       :config:`debris.transport.terminus_removal`
     - Switch the gravitational transport and the removal at the margin on or off.
   * - :config:`debris.transport.input.file`, :config:`debris.transport.input.periodic`
     - The debris input forcing (leave the file name empty for no input).
   * - :config:`debris.transport.cover_fraction_coefficient`
     - Coefficient `C` of the debris-covered area fraction `1 - \exp(-C h_d)`.

The model state (:var:`debris_thickness` and the englacial concentration
:var:`englacial_debris_concentration`, kg m\ :sup:`-3`) is saved to output files and read
when re-starting; when bootstrapping both are read from the input file if present and set
to zero otherwise.

Diagnostics: :var:`debris_thickness`, :var:`englacial_debris_concentration`,
:var:`englacial_debris_column_mass`, :var:`debris_cover_fraction`,
:var:`debris_input_rate`, :var:`debris_melt_out_rate`, :var:`debris_removal_rate`,
:var:`debris_surface_velocity`, :var:`debris_gravitational_flux`; scalar time series
:var:`englacial_debris_mass`, :var:`supraglacial_debris_mass`, :var:`total_debris_mass`,
:var:`debris_input_mass_flux`, :var:`debris_melt_out_mass_flux`,
:var:`debris_output_mass_flux`, :var:`debris_lost_mass_flux` (debris removed with basal
melt or in cells that became ice-free) and :var:`debris_mass_conservation_error` (the
debris mass minus its initial value and the net of inputs, outputs and losses since the
start of the run).

See ``examples/debris/`` for the idealized valley glacier of the paper.

.. _sec-debris-melt-enhancement:

Effect of the debris cover on melt
++++++++++++++++++++++++++++++++++

:|options|: :config:`debris.ice_melt_enhancement.model`
:|variables|: :var:`debris_melt_factor` (with the ``given`` model)
:|implementation|: ``pism::debris::IceMeltEnhancement``

The dimensionless factor :var:`ice_melt_enhancement` multiplying the clean-ice melt rate
is computed from the debris thickness using the Østrem curve of
:cite:`VerhaegenHuybrechts2026` (equation 14; model ``verhaegen``), read from
:config:`debris.ice_melt_enhancement.file` (model ``given``), or set to one (model
``none``).
