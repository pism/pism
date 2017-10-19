.. include:: ../../../global.txt

.. _sec-energy:

Modeling conservation of energy
-------------------------------

In normal use PISM solves the conservation of energy problem within the ice, the thin
subglacial layer, and a layer of thermal bedrock. For the ice and the subglacial layer it
uses an enthalpy-based scheme :cite:`AschwandenBuelerKhroulevBlatter` which allows the energy
to be conserved even when the temperature is at the pressure-melting point.

Ice at the melting point is called "temperate" ice. Part of the thermal energy of
temperate ice is in the latent heat of the liquid water stored between the crystals of the
temperate ice. Part of the thermal energy of the whole glacier is in the latent heat of
the liquid water under the glacier. The enthalpy scheme correctly models these storehouses
of thermal energy, and thus it allows polythermal and fully-temperate glaciers to be
modeled :cite:`AschwandenBlatter`.

The state of the full conservation of energy model includes the 3D ``enthalpy`` variable
plus the 2D ``bwat`` and ``tillwat`` subglacial hydrology state variables (subsection
:ref:`sec-subhydro`), all of which are seen in output files. The important basal melt rate
computation involves all of these energy state variables, because the basal melt rate
(``bmelt`` in output files) comes from conserving energy across the ice-bedrock layer
:cite:`AschwandenBuelerKhroulevBlatter`. Fields ``temp``, ``liqfrac``, and ``temp_pa`` seen in
output files are all actually diagnostic outputs because all of these can be recovered
from the enthalpy and the ice geometry.

Because this part of PISM is just a conservation law, there is little need for the user to
worry about controlling it. If desired, however, conservation of energy can be turned off
entirely with :opt:`-energy none`. The default enthalpy-based conservation of energy model
(i.e. ``-energy enthalpy``) can be replaced by the temperature-based (i.e. "cold ice")
method used in :cite:`BBssasliding` and verified in :cite:`BBL` by setting option :opt:`-energy
cold`.

The thermal bedrock layer model is turned off by setting ``-Mbz 1`` (i.e. zero spaces)
while it is turned on by choosing a depth and number of points, as in ``-Lbz 1000 -Mbz
21``, for example, which gives a layer depth of 1000 m and grid spaces of 50 m (=
1000/20). The input geothermal flux (``bheatflx`` in output files) is applied at the
bottom of the bedrock thermal layer if such a layer is present and otherwise it is applied
at the base of the ice.
