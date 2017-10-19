.. include:: ../../../global.txt

.. _sec-gradient:

Surface gradient method
-----------------------

PISM computes surface gradients to determine the "driving stress"

.. math::

   (\tau_{d,x},\tau_{d,y}) = - \rho g H \nabla h,

where `H` is the ice thickness, and `h` is the ice surface elevation. The driving stress
enters into both the SIA and SSA stress balances, but in the former the driving stress is
needed on a staggered grid, while in the latter the driving stress is needed on the
regular grid.

Surface gradients are computed by finite differences in several slightly-different ways.
There are options for choosing which method to use **in the SIA model**, but to the best
of our knowledge there is no theoretical advice on the best, most robust mechanism.

The SSA model uses centered finite differences, switching from centered to one-sided near
the ice margin and does not recognize choices other than :opt:`-gradient eta`.

There are three :opt:`-gradient` methods in PISM:

.. list-table:: Options controlling the surface gradient computation in the SIA code
   :name: tab-sia-gradient
   :header-rows: 1
   :widths: 1,3

   * - Option
     - Description

   * - :opt:`-gradient mahaffy`
     - This most "standard" way computes the surface slope onto the staggered grid for the
       SIA :cite:`Mahaffy`. It makes `O(\Delta x^2,\Delta y^2)` errors. For computations
       of driving stress on the regular grid, centered differencing is used instead.

   * - :opt:`-gradient haseloff`
     - This is the default method. It only differs from ``mahaffy`` at ice-margin
       locations, where the slope is approximated using one-sided finite differences in
       cases where an adjacent ice-free bedrock surface elevation is above the ice
       elevation.

   * - :opt:`-gradient eta`
     - In this method we first transform the thickness `H` by `\eta = H^{(2n+2)/n}` and
       then differentiate the sum of the thickness and the bed using centered differences:

       .. math::

          \nabla h = \nabla H + \nabla b = \frac{n}{(2n+2)}
          \eta^{(-n-2)/(2n+2)} \nabla \eta + \nabla b.

       Here `b` is the bed elevation and `h` is the surface elevation. This transformation
       sometimes has the benefits that the surface values of the horizontal velocity and
       vertical velocity, and the driving stress, are better behaved near the margin. See
       :cite:`BLKCB` for technical explanation of this transformation and compare
       :cite:`SaitoMargin`. The actual finite difference schemes applied to compute the
       surface slope are similar to option ``mahaffy``.

.. note:: The :opt:`-gradient eta` may improve the model performance near *grounded*
          margins but should not be used in simulations of marine ice sheets.
