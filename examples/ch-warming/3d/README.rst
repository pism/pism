A minimal 3D setup using the parameterization of cryo-hydrological warming
==========================================================================

Run `run.sh` to tests a minimal cryo-hydrological warming setup.

This example repeats the one-column test setup using `pismr -regional` instead of PISM's
Python bindings.

.. note::

   The option `-regional` *is* required at this point: the cryo-hydrologic warming
   parameterization is implemented as a part of PISM's regional code.

Set `energy.ch_warming.enabled` to enable this parameterization.

To control it, adjust

- `energy.ch_warming.average_channel_spacing` (`R` in the paper)
- `energy.ch_warming.residual_water_fraction` (residual water fraction in the
  cryo-hydrological system during (and at the end) of the melt season)
- `energy.ch_warming.temperate_ice_thermal_conductivity_ratio` (the ratio of thermal
  conductivities of temperate and cold ice in the cryo-hydrological system). This
  parameter is set to `1` to make temperate ice more conductive (as in the paper) compared
  to PISM's model.
