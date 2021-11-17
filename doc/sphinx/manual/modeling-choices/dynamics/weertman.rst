.. include:: ../../../global.txt

.. _sec-weertman:

Weertman-style sliding law
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning::

   This kind of sliding is, in general, a bad idea. We implement it to simplify
   comparisons of the "hybrid" model mentioned above to older studies using this
   parameterization.

The "Weertman-type sliding law" (:cite:`GreveBlatter2009`, equations 5.35 and 5.91) has
the form

.. math::

   \mathbf{u}_s =
   \begin{cases}
   \mathbf{0}, & T_b < T_m, \\
   -C_b(\rho g H)^{p-q}|\nabla h|^{p-1}\nabla h, & T_b = T_m,
   \end{cases}

`T_b` is the ice temperature, and `T_m` is the pressure-melting temperature. The constant
`C_b` and exponents `p` and `q` are tuning parameters.

The particular form implemented in PISM comes from equation 5 in :cite:`Tomkin2007`:

.. math::
   :label: eq-weertman-sliding

   \mathbf{u}_s = -\frac{2 A_s \beta_c (\rho g H)^{n}}{N - P} |\nabla h|^{n-1} \nabla h.

.. list-table:: Notation used in :eq:`eq-weertman-sliding`
   :name: tab-weertman-notation
   :header-rows: 1
   :widths: 1,9

   * - Variable
     - Meaning

   * - `H`
     - ice thickness

   * - `h`
     - ice surface elevation

   * - `n`
     - flow law exponent

   * - `g`
     - acceleration due to gravity

   * - `\rho`
     - ice density

   * - `N`
     - ice overburden pressure, `N = \rho g H`

   * - `P`
     - basal water pressure

   * - `A_s`
     - sliding parameter

   * - `\beta_c`
     - "constriction parameter" capturing the effect of valley walls on the flow;
       set to `1` in this implementation

We assume that the basal water pressure is a given constant fraction of the overburden
pressure: `P = k N`. This simplifies :eq:`eq-weertman-sliding` to

.. math::

   \mathbf{u}_s = -\frac{2 A_s}{1 - k} ( \rho g H\, |\nabla h| )^{n-1} \nabla h.

This parameterization is used for grounded ice *where the base of the ice is temperate*.

To enable, use :opt:`-stress_balance weertman_sliding` (this results in constant-in-depth
ice velocity) or :opt:`-stress_balance weertman_sliding+sia` to use this parameterization
as a sliding law with the deformational flow modeled using the SIA model.

Use configuration parameters :config:`stress_balance.weertman_sliding.k` and
:config:`stress_balance.weertman_sliding.A` to set `k` and `A_s`, respectively. Default
values come from :cite:`Tomkin2007`.

Parameters
##########

Prefix: ``stress_balance.weertman_sliding.``

.. pism-parameters::
   :prefix: stress_balance.weertman_sliding.
