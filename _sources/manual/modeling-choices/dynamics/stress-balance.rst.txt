.. include:: ../../../global.txt

.. _sec-stressbalance:

Choosing the stress balance
---------------------------

The basic stress balance used for all grounded ice in PISM is the non-sliding,
thermomechanically-coupled SIA :cite:`BBL`. For the vast majority of most ice sheets, as
measured by area or volume, this is an appropriate model, which is an `O(\epsilon^2)`
approximation to the Stokes model if `\epsilon` is the depth-to-length ratio of the ice
sheet :cite:`Fowler`.

The shallow shelf approximation (SSA) stress balance applies to floating ice. See the Ross
ice shelf example in section :ref:`sec-ross` for an example in which the SSA is only
applied to floating ice.

In PISM the SSA is also used to describe the sliding of grounded ice and the formation of
ice streams :cite:`BBssasliding`. Specifically for the SSA with "plastic" (Coulomb
friction) basal resistance, the locations of ice streams are determined as part of a free
boundary problem of Schoof :cite:`SchoofStream`, a model for emergent ice streams within a
ice sheet and ice shelf system. This model explains ice streams through a combination of
plastic till failure and SSA stress balance.

This SSA description of ice streams is the preferred "sliding law" for the SIA
:cite:`BBssasliding`, :cite:`Winkelmannetal2011`. The SSA should be combined with the SIA,
in this way, in preference to classical SIA sliding laws which make the sliding velocity
of ice a local function of the basal value of the driving stress. The resulting
combination of SIA and SSA is a "hybrid" approximation of the Stokes model
:cite:`Winkelmannetal2011`. Option ``-stress_balance ssa+sia`` turns on this "hybrid"
model. In this use of the SSA as a sliding law, floating ice is also subject to the SSA.

In addition to this, PISM includes an implementation of the first order approximation of
Stokes equations due to Blatter (``-stress_balance blatter``, :cite:`Blatter`,
:cite:`Pattyn03`).

All stress balance options *except* for the first order approximation correspond to two
basic choices:

- modeling basal sliding, and
- modeling of ice velocity within an ice column.

PISM supports the following stress balance choices, controlled using
:config:`stress_balance.model` (option :opt:`-stress_balance`):

#. ``none``: no sliding, ice velocity is constant in each column. This
   equivalent to disabling ice flow completely.

#. ``prescribed_sliding``: Use the constant-in-time prescribed sliding velocity field read
   from a file set using :config:`stress_balance.prescribed_sliding.file`, variables
   ``ubar`` and ``vbar``. Horizontal ice velocity is constant throughout ice columns.

#. ``ssa``: Use the :ref:`sec-ssa` model exclusively. Horizontal ice
   velocity is constant throughout ice columns.

#. ``weertman_sliding``: basal sliding is approximated using the
   :ref:`sec-weertman`, ice velocity is constant throughout ice columns.

#. ``sia`` (*default*): no sliding; ice velocity within the column is approximated using
   the :ref:`sec-sia`. Floating ice does not flow, so this model is not recommended for
   marine ice sheets.

#. ``prescribed_sliding+sia``: basal ice velocity is read from an input
   file and held constant, ice velocity within the column is approximated using the
   :ref:`sec-sia`.

#. ``ssa+sia``: use :ref:`sec-ssa` as a sliding law with a plastic or
   pseudo-plastic till, combining it with the :ref:`sec-sia` according to the combination
   in :cite:`Winkelmannetal2011`; similar to :cite:`BBssasliding`. Floating ice uses SSA
   only. *This "hybrid" stress balance is the recommended sliding law for the SIA.*

#. ``weertman_sliding+sia``: basal sliding is approximated using the
   :ref:`sec-weertman`, ice velocity within the column is approximated using the
   :ref:`sec-sia`.

#. ``blatter``: use :ref:`sec-blatter`.

Please see the following sections for details.

.. toctree::

   ssa.rst

   sia.rst

   weertman.rst

   blatter.rst
