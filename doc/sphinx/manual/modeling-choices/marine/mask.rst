.. include:: ../../../global.txt

.. _sec-floatmask:

Flotation criterion, mask, and sea level
----------------------------------------

The most basic decision about marine ice sheet dynamics made internally by PISM is whether
a ice-filled grid cell is floating. That is, PISM applies the "flotation criterion"
:cite:`Winkelmannetal2011` at every time step and at every grid location to determine whether
the ice is floating on the ocean or not. The result is stored in the ``mask`` variable.
The ``mask`` variable has ``pism_intent`` = ``diagnostic``, and thus it does *not* need to
be included in the input file set using the ``-i`` option.

The possible values of the ``mask`` are given in :numref:`tab-maskvals`. The mask
does not *by itself* determine ice dynamics. For instance, even when ice is floating (mask
value ``MASK_FLOATING``), the user must turn on the usual choice for ice shelf dynamics,
namely the SSA stress balance, by using options :opt:`-stress_balance ssa` or
:opt:`-stress_balance ssa+sia`.

.. list-table:: The PISM mask, in combination with user options, determines the dynamical
                model.
   :name: tab-maskvals
   :header-rows: 1
   :widths: 1,1

   * - Mask value
     - Meaning

   * - 0 = ``MASK_ICE_FREE_BEDROCK``
     - ice free bedrock 

   * - 2 = ``MASK_GROUNDED``
     - ice is grounded 

   * - 3 = ``MASK_FLOATING``
     - ice is floating (the SIA is never applied; the SSA is applied if the ``ssa`` or
       ``ssa+sia`` stress balance model is selected

   * - 4 = ``MASK_ICE_FREE_OCEAN``
     - ice-free ocean 

Assuming that the geometry of the ice is allowed to evolve (which can be turned off by
option ``-no_mass``), and assuming an ocean exists so that a sea level is used in the
flotation criterion (which can be turned off by option :opt:`-dry`), then at each time
step the mask will be updated.
