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

- setting :config:`geometry.update.enabled` to "false" disables the mass-continuity
  (conservation of mass) step
- setting :config:`energy.model` to "none" disables the conservation of energy computation
- setting :config:`energy.model` to "cold" makes PISM use temperature instead of enthalpy
  in the energy conservation code
- setting :config:`stress_balance.model` to "none" disables the stress balance computation
  (useful for testing surface mass balance inputs)
