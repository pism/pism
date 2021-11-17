.. include:: ../../../global.txt

.. _sec-sia:

Shallow ice approximation (SIA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

   The explicit time stepping of the mass continuity equation in the case of the SIA flow
   comes with a severe restriction on time step length:

   .. math::
      :label: eq-sia-max-dt

      \dt \le \frac{2 R}{D\left( 1/\dx^2 + 1/\dy^2 \right)}

   Here `D` is the maximum diffusivity of the SIA flow and `R` is
   :config:`time_stepping.adaptive_ratio`, a tuning parameter that further reduces the
   maximum allowed time step length.

   The maximum diffusivity `D` may be achieved at an isolated grid point near the ice
   margin. In this case it might make sense to limit the diffusivity of the SIA flow,
   sacrificing accuracy at a few grid points to increase time step length and reduce the
   computational cost. Set :config:`stress_balance.sia.limit_diffusivity` to enable this
   mechanism.

   When :config:`stress_balance.sia.limit_diffusivity` is ``false`` PISM stops as soon as
   the SIA diffusivity at any grid point exceeds
   :config:`stress_balance.sia.max_diffusivity`. We do this to make it easier to detect
   problematic model configurations: in many cases it does not make sense to continue a
   simulation if `D` is very large.

.. _sec-sia-gradient:

Surface gradient method
-----------------------

PISM computes surface gradients to determine the "driving stress"

.. math::

   (\tau_{d,x},\tau_{d,y}) = - \rho g H \nabla h,

where `H` is the ice thickness, and `h` is the ice surface elevation.

In the SIA model surface gradients at *staggered* grid locations are computed using one of
the following three finite-difference approximations (selected using
:config:`stress_balance.sia.surface_gradient_method`):

#. ``mahaffy``: This most "standard" way computes the surface slope onto the staggered
   grid for the SIA :cite:`Mahaffy`. It makes `O(\dx^2,\dy^2)` errors.

#. ``haseloff``: This is the default method. It only differs from ``mahaffy`` at
   ice-margin locations, where the slope is approximated using one-sided finite
   differences in cases where an adjacent ice-free bedrock surface elevation is above the
   ice elevation.

#. ``eta``: In this method we first transform the thickness `H` by `\eta = H^{(2n+2)/n}`
   and then differentiate the sum of the thickness and the bed using centered differences:

   .. math::

      \nabla h = \nabla H + \nabla b = \frac{n}{(2n+2)}
      \eta^{(-n-2)/(2n+2)} \nabla \eta + \nabla b.

   Here `b` is the bed elevation, `h` is the surface elevation, and `n` is the Glen
   exponent. This transformation sometimes has the benefits that the surface values of the
   horizontal velocity and vertical velocity, and the driving stress, are better behaved
   near the margin. See :cite:`BLKCB` for technical explanation of this transformation and
   compare :cite:`SaitoMargin`. The actual finite difference schemes applied to compute
   the surface slope are similar to option ``mahaffy``.

   .. note:: This method may improve the model performance near *grounded* margins but
             should not be used in simulations of marine ice sheets.

To the best of our knowledge there is no theoretical advice on the best, most robust
mechanism.

.. _sec-bedsmooth:

Parameterization of bed roughness
---------------------------------

Schoof :cite:`Schoofbasaltopg2003` describes how to alter the SIA stress balance to model
ice flow over bumpy bedrock topgraphy. One computes the amount by which bumpy topography
lowers the SIA diffusivity. An internal quantity used in this method is a smoothed version
of the bedrock topography. As a practical matter for PISM, this theory improves the SIA's
ability to handle bed roughness because it parameterizes the effects of "higher-order"
stresses which act on the ice as it flows over bumps. For additional technical description
of PISM's implementation, see :ref:`sec-bed-roughness`.

There are two associated parameters:

- :config:`stress_balance.sia.bed_smoother.range` gives the half-width of the square
  smoothing domain in meters. If zero is given then the mechanism is turned off. The
  mechanism is on by default using executable ``pismr``, with the half-width set to 5 km,
  giving Schoof's recommended smoothing size of 10 km :cite:`Schoofbasaltopg2003`.
- :config:`stress_balance.sia.bed_smoother.theta_min` is the minimum value of `\theta` in
  the parameterization.

This mechanism is turned off by default in ``pismv``.

Under the default :config:`output.size` (``medium``), PISM writes fields :var:`topgsmooth`
and :var:`schoofs_theta` from this mechanism. The thickness relative to the smoothed
bedrock elevation, namely :var:`topgsmooth`, is the difference between the unsmoothed
surface elevation and the smoothed bedrock elevation. It is *only used internally by this
mechanism*, to compute a modified value of the diffusivity; the rest of PISM does not use
this or any other smoothed bed. The field :var:`schoofs_theta` is a number `\theta`
between :config:`stress_balance.sia.bed_smoother.theta_min` (usually zero) and `1`, with
values significantly below one indicating a reduction in diffusivity, essentially a drag
coefficient, from bumpy bed.

.. _sec-sia-age-coupling:

Coupling to the age of the ice
------------------------------

The age of the ice can be used in two parameterizations in the SIA stress balance model:

#. Ice grain size parameterization based on data from :cite:`DeLaChapelleEtAl98` and
   :cite:`LipenkovEtAl89` (Vostok core data). In PISM, only the Goldsby-Kohlstedt flow law
   (see :ref:`sec-rheology`) uses the grain size.

   Set :config:`stress_balance.sia.grain_size_age_coupling` to enable this parameterization.

#. The flow enhancement factor can be coupled to the age of the ice as in
   :cite:`Greve97Greenland`: `e` in :eq:`eq-sia-enhancement` is set to

   - :config:`stress_balance.sia.enhancement_factor_interglacial` during Eemian and Holocene,
   - :config:`stress_balance.sia.enhancement_factor` otherwise.

   See :config:`time.eemian_start`, :config:`time.eemian_end`, and
   :config:`time.holocene_start`.

   Set :config:`stress_balance.sia.e_age_coupling` to enable this parameterization.

.. _sec-sia-parameters:

Parameters
----------

Prefix: ``stress_balance.sia.``

.. pism-parameters::
   :prefix: stress_balance.sia.
