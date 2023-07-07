.. include:: ../global.txt

.. _sec-cf-standard-names:

CF standard names used by PISM
==============================

Existing standard names
-----------------------

We start by listing standard names from the CF Standard Name Table. The subset here is a
small subset of the `table <CF-standard-names_>`_; we list only

- those with "land_ice" in the name and
- those currently used by PISM

The existing names starting with "land_ice" are believed to have all been submitted by
Magnus Hagdorn to the CF committee circa 2003. The `SeaRISE assessment process <searise_>`_
now has a `wiki on CF standard name use <cf-names-for-glaciology_>`_, which to a
significant extent duplicates content regarding proposed names on this page. That wiki is
an evolving community standard, and it supercedes this page when it comes to actual
evolving standards.

Go to the `CF Conventions <cf-names-proposed_>`_ page for the list of proposed standard
names under consideration.

Because of the use of UDUNITS_, PISM input files do not have to have fields already in the
canonical units. Rather, the units attribute has to be valid for UDUNITS conversion into
the canonical units. Generally within PISM, the canonical units are used internally.

.. list-table:: CF standard names used by PISM
   :name: tab-standard-names
   :header-rows: 1

   * - CF standard name
     - Canonical units (SI)

   * - ``bedrock_altitude``
     - m
   * - ``floating_ice_sheet_area_fraction``
     - 1
   * - ``grounded_ice_sheet_area_fraction``
     - 1
   * - ``land_ice_area_fraction``
     - 1
   * - ``land_ice_basal_melt_rate``
     - m s-1
   * - ``land_ice_basal_temperature``
     - K
   * - ``land_ice_basal_upward_velocity``
     - m s-1
   * - ``land_ice_basal_x_velocity``
     - m s-1
   * - ``land_ice_basal_y_velocity``
     - m s-1
   * - ``land_ice_calving_rate``
     - m s-1
   * - ``land_ice_lwe_basal_melt_rate``
     - m s-1
   * - ``land_ice_lwe_calving_rate``
     - m s-1
   * - ``land_ice_lwe_surface_specific_mass_balance``
     - m s-1
   * - ``land_ice_sigma_coordinate``
     - 1
   * - ``land_ice_specific_mass_flux_due_to_calving_and_ice_front_melting``
     - kg m-2 s-1
   * - ``land_ice_surface_specific_mass_balance_flux``
     - kg m-2 s-1
   * - ``land_ice_surface_upward_velocity``
     - m s-1
   * - ``land_ice_surface_x_velocity``
     - m s-1
   * - ``land_ice_surface_y_velocity``
     - m s-1
   * - ``land_ice_temperature``
     - K
   * - ``land_ice_thickness``
     - m
   * - ``land_ice_vertical_mean_x_velocity``
     - m s-1
   * - ``land_ice_vertical_mean_y_velocity``
     - m s-1
   * - ``land_ice_x_velocity``
     - m s-1
   * - ``land_ice_y_velocity``
     - m s-1
   * - ``latitude``
     - degree_north
   * - ``longitude``
     - degree_east
   * - ``magnitude_of_land_ice_basal_drag``
     - Pa
   * - ``projection_x_coordinate``
     - m
   * - ``projection_y_coordinate``
     - m
   * - ``surface_altitude``
     - m
   * - ``temperature_at_ground_level_in_snow_or_firn``
     - K
   * - ``tendency_of_bedrock_altitude``
     - m s-1
   * - ``tendency_of_land_ice_mass_due_to_basal_mass_balance``
     - kg s-1
   * - ``tendency_of_land_ice_thickness``
     - m s-1
   * - ``upward_geothermal_heat_flux_at_ground_level``
     - W m-2


*Proposed* standard names
-------------------------

These are *unofficially* proposed by Bueler and Aschwanden, for now.

.. list-table:: Desired CF standard names
   :header-rows: 1
   :widths: 7,2,6

   * - Proposed name
     - Canonical units (SI)
     - Comments

   * - ``ice_shelf_basal_specific_mass_balance``
     - m s-1
     - positive is loss of ice shelf mass (i.e. use outward normal from ice shelf)
   * - ``ice_shelf_basal_temperature``
     - K
     - absolute (not pressure-adjusted) temperature
   * - ``land_ice_age``
     - s
     -
   * - ``land_ice_basal_frictional_heating``
     - W m-2
     -
   * - ``land_ice_basal_material_yield_stress``
     - Pa
     -
   * - ``land_ice_basal_material_friction_angle``
     - degree
     - majority of standard names with "angle" use canonical units "degree"
   * - ``land_ice_surface_temperature_below_firn``
     - K
     -
   * - ``land_ice_upward_velocity``
     - m s-1
     - compare to CF names "upward_air_velocity" and "upward_sea_water_velocity"
   * - ``lithosphere_temperature``
     - K
     -
   * - ``upward_geothermal_flux_in_lithosphere``
     - W m-2
     - typically applied at depth in lithosphere; compare to
       "upward_geothermal_heat_flux_at_ground_level"
   * - ``land_ice_specific_enthalpy``
     - J kg-1
     - enthalpy is defined in PISM to be sensible plus latent heat, plus potential energy
       of pressure; there is a nontrivial issue of the scaling; the enthalpy value for
       273.15 K (cold) ice at atmospheric pressure is a possible standard
   * - ``land_ice_liquid_fraction``
     - 1
     - liquid water fraction in ice, a pure number between 0 and 1; a diagnostic function
       of enthalpy


Final technical notes
---------------------

- PISM also uses attributes ``grid_mapping = "mapping" ;`` and ``coordinates = "lat lon";``
  on output variables that depend on ``y,x``.
- Because PISM uses UDUNITS, it will write some variables in "human-friendly" units
  instead of the SI units listed above, for instance velocities in ``m year-1`` instead of
  ``m s-1``. This is allowed under CF. When PISM reads such a field from a NetCDF file,
  the conversion is handled automatically by UDUNITS_.
