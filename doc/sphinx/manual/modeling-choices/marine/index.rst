.. include:: ../../../global.txt

.. _sec-marine:

Marine ice sheet modeling
=========================

PISM is often used to model whole ice sheets surrounded by ocean, with attached floating
ice shelves, or smaller regions like outlet glaciers flowing into embayments and possibly
generating floating tongues. This section explains the geometry and stress balance
mechanisms in PISM that apply to floating ice, at the vertical calving faces of floating
ice, or at marine grounding lines. The physics at calving fronts is very different from
elsewhere on an ice sheet, because the flow is nothing like the lubrication flow addressed
by the SIA, and nor is the physics like the sliding flow in the interior of an ice domain.
The needed physics at the calving front can be thought of as boundary condition
modifications to the mass continuity equation and to the SSA stress balance equation. The
physics of grounding lines are substantially handled by recovering sub-grid information
through interpolation.

.. toctree::

   pik.rst

   mask.rst

   calving.rst

   melange.rst
