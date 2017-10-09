.. include:: ../../global.txt

.. _sec-turning-off:

Disabling sub-models
--------------------

Certain major model components, unlike more peripheral ones like bed deformation or
calving, are "on" by default. They do not need to be turned on explicitly. For example,
the SIA computation is so common that it would be a hassle to require an option to turn it
on every time you need it.

But sometimes one wants to disable particular components, during model spin-up, for
example. PISM has the following "off" switches:

- :opt:`-no_mass` disables the mass-continuity (conservation of mass) step
- :opt:`-energy none` disables the conservation of energy computation
- :opt:`-energy cold` makes PISM use temperature instead of enthalpy in the energy
  conservation code
- :opt:`-stress_balance none` disables the stress balance computation (useful for testing
  surface mass balance inputs)
- :opt:`-dry` essentially disables ocean models: ice is always considered to be grounded,
  the sub-shelf melt rate and temperature is not used, and the calving-front boundary
  condition is computed ignoring the water pressure exerted on the vertical face at a
  (possibly submerged) terminus.
