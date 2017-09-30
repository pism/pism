
.. DO NOT EDIT. This file was generated using list_diagnostics.py.

.. _sec-diagnostics-list:

List of PISM's diagnostics
==========================

.. _sec-extra_vars:

Spatially-variable fields
-------------------------

#. ``basal_mass_flux_floating``

   :Units: kg m-2 year-1
   :Description: average basal mass flux over the reporting interval (floating areas)
   :Standard name: ---
   :Comment: positive flux corresponds to ice gain

#. ``basal_mass_flux_grounded``

   :Units: kg m-2 year-1
   :Description: average basal mass flux over the reporting interval (grounded areas)
   :Standard name: ---
   :Comment: positive flux corresponds to ice gain

#. ``basal_melt_rate_grounded``

   :Units: m year-1
   :Description: ice basal melt rate from energy conservation, in ice thickness per time (valid in grounded areas)
   :Standard name: ---
   :Comment: positive basal melt rate corresponds to ice loss

#. ``bedtoptemp``

   :Units: Kelvin
   :Description: temperature at the top surface of the bedrock thermal layer
   :Standard name: ---

#. ``beta``

   :Units: Pa s / m
   :Description: basal drag coefficient
   :Standard name: ---

#. ``bfrict``

   :Units: W m-2
   :Description: basal frictional heating
   :Standard name: ---

#. ``bheatflx``

   :Units: mW m-2
   :Description: upward geothermal flux at the bottom bedrock surface
   :Standard name: ---
   :Comment: positive values correspond to an upward flux

#. ``bmelt``

   :Units: m year-1
   :Description: ice basal melt rate from energy conservation and subshelf melt, in ice thickness per time
   :Standard name: ``land_ice_basal_melt_rate``
   :Comment: positive basal melt rate corresponds to ice loss

#. ``bwat``

   :Units: m
   :Description: thickness of transportable water in subglacial layer
   :Standard name: ---

#. ``bwp``

   :Units: Pa
   :Description: pressure of transportable water in subglacial layer
   :Standard name: ---

#. ``bwprel``

   :Units: ---
   :Description: pressure of transportable water in subglacial layer as fraction of the overburden pressure
   :Standard name: ---

#. ``cell_area``

   :Units: km2
   :Description: cell areas
   :Standard name: ---
   :Comment: values are equal to dx*dy if projection parameters are not available; otherwise WGS84 ellipsoid is used

#. ``cell_grounded_fraction``

   :Units: ---
   :Description: fractional grounded/floating mask (floating=0, grounded=1)
   :Standard name: ---

#. ``climatic_mass_balance``

   :Units: kg m-2 year-1
   :Description: surface mass balance (accumulation/ablation) rate
   :Standard name: ``land_ice_surface_specific_mass_balance_flux``

#. ``cts``

   :Units: ---
   :Description: cts = E/E_s(p), so cold-temperate transition surface is at cts = 1
   :Standard name: ---

#. ``dHdt``

   :Units: m year-1
   :Description: ice thickness rate of change
   :Standard name: ``tendency_of_land_ice_thickness``

#. ``dbdt``

   :Units: mm year-1
   :Description: bedrock uplift rate
   :Standard name: ``tendency_of_bedrock_altitude``

#. ``deviatoric_stresses``

   - ``sigma_xx``

     :Units: Pa
     :Description: deviatoric stress in X direction
     :Standard name: ---

   - ``sigma_yy``

     :Units: Pa
     :Description: deviatoric stress in Y direction
     :Standard name: ---

   - ``sigma_xy``

     :Units: Pa
     :Description: deviatoric shear stress
     :Standard name: ---

#. ``diffusivity``

   :Units: m2 s-1
   :Description: diffusivity of SIA mass continuity equation
   :Standard name: ---

#. ``diffusivity_staggered``

   - ``diffusivity_i``

     :Units: m2 s-1
     :Description: diffusivity of SIA mass continuity equation on the staggered grid (i-offset)
     :Standard name: ---

   - ``diffusivity_j``

     :Units: m2 s-1
     :Description: diffusivity of SIA mass continuity equation on the staggered grid (j-offset)
     :Standard name: ---

#. ``effbwp``

   :Units: Pa
   :Description: effective pressure of transportable water in subglacial layer (overburden pressure minus water pressure)
   :Standard name: ---

#. ``effective_viscosity``

   :Units: kPascal second
   :Description: effective viscosity of ice
   :Standard name: ---

#. ``enthalpy``

   :Units: J kg-1
   :Description: ice enthalpy (includes sensible heat, latent heat, pressure)
   :Standard name: ---

#. ``enthalpybase``

   :Units: J kg-1
   :Description: ice enthalpy at the base of ice
   :Standard name: ---

#. ``enthalpysurf``

   :Units: J kg-1
   :Description: ice enthalpy at 1m below the ice surface
   :Standard name: ---

#. ``flux``

   - ``uflux``

     :Units: m2 year-1
     :Description: Vertically integrated horizontal flux of ice in the X direction
     :Standard name: ---

   - ``vflux``

     :Units: m2 year-1
     :Description: Vertically integrated horizontal flux of ice in the Y direction
     :Standard name: ---

#. ``flux_divergence``

   :Units: m year-1
   :Description: flux divergence
   :Standard name: ---

#. ``flux_mag``

   :Units: m2 year-1
   :Description: magnitude of vertically-integrated horizontal flux of ice
   :Standard name: ---

#. ``flux_staggered``

   :Units: m2 year-1
   :Description: fluxes through cell interfaces (sides) on the staggered grid
   :Standard name: ---

#. ``h_x``

   - ``h_x_i``

     :Units: ---
     :Description: the x-component of the surface gradient, i-offset
     :Standard name: ---

   - ``h_x_j``

     :Units: ---
     :Description: the x-component of the surface gradient, j-offset
     :Standard name: ---

#. ``h_y``

   - ``h_y_i``

     :Units: ---
     :Description: the y-component of the surface gradient, i-offset
     :Standard name: ---

   - ``h_y_j``

     :Units: ---
     :Description: the y-component of the surface gradient, j-offset
     :Standard name: ---

#. ``hardav``

   :Units: Pa s0.333333
   :Description: vertical average of ice hardness
   :Standard name: ---

#. ``hardness``

   :Units: Pa s0.333333
   :Description: ice hardness computed using the SIA flow law
   :Standard name: ---

#. ``height_above_flotation``

   :Units: m
   :Description: ice thickness in excess of the maximum floating ice thickness
   :Standard name: ---
   :Comment: shows how close to floatation the ice is at a given location

#. ``hfgeoubed``

   :Units: mW m-2
   :Description: upward geothermal flux at the top bedrock surface
   :Standard name: ``upward_geothermal_heat_flux_at_ground_level``
   :Comment: positive values correspond to an upward flux

#. ``hydrobmelt``

   :Units: m year-1
   :Description: the version of bmelt seen by the hydrology model
   :Standard name: ---

#. ``hydroinput``

   :Units: m year-1
   :Description: total water input into subglacial hydrology layer
   :Standard name: ---

#. ``ice_area_specific_volume``

   :Units: m3/m2
   :Description: ice-volume-per-area in partially-filled grid cells
   :Standard name: ---
   :Comment: this variable represents the amount of ice in a partially-filled cell and not the corresponding geometry, so thinking about it as 'thickness' is not helpful

#. ``ice_mass``

   :Units: kg
   :Description: mass per cell
   :Standard name: ---

#. ``ice_surface_liquid_water_fraction``

   :Units: 1
   :Description: ice liquid water fraction at the ice surface
   :Standard name: ---

#. ``ice_surface_temp``

   :Units: Kelvin
   :Description: ice temperature at the ice surface
   :Standard name: ---

#. ``lat``

   :Units: degree_north
   :Description: latitude
   :Standard name: ``latitude``

#. ``liqfrac``

   :Units: 1
   :Description: liquid water fraction in ice (between 0 and 1)
   :Standard name: ---

#. ``lon``

   :Units: degree_east
   :Description: longitude
   :Standard name: ``longitude``

#. ``mask``

   :Units: ---
   :Description: ice-type (ice-free/grounded/floating/ocean) integer mask
   :Standard name: ---

#. ``melange_back_pressure_fraction``

   :Units: 1
   :Description: dimensionless pressure fraction at calving fronts due to presence of melange 
   :Standard name: ---

#. ``ocean_pressure_difference``

   :Units: ---
   :Description: ocean pressure difference at calving fronts
   :Standard name: ---

#. ``pressure``

   :Units: Pa
   :Description: pressure in ice (hydrostatic)
   :Standard name: ---

#. ``rank``

   :Units: 1
   :Description: processor rank
   :Standard name: ---

#. ``schoofs_theta``

   :Units: 1
   :Description: multiplier 'theta' in Schoof's (2003) theory of bed roughness in SIA
   :Standard name: ---

#. ``sea_level``

   :Units: meters
   :Description: sea level elevation, relative to the geoid
   :Standard name: ---

#. ``sftflf``

   :Units: 1
   :Description: fraction of a grid cell covered by floating ice
   :Standard name: ``floating_ice_sheet_area_fraction``

#. ``sftgif``

   :Units: 1
   :Description: fraction of a grid cell covered by ice (grounded or floating)
   :Standard name: ``land_ice_area_fraction``

#. ``sftgrf``

   :Units: 1
   :Description: fraction of a grid cell covered by grounded ice
   :Standard name: ``grounded_ice_sheet_area_fraction``

#. ``shelfbmassflux``

   :Units: kg m-2 s-1
   :Description: mass flux at the basal surface of ice shelves
   :Standard name: ---

#. ``shelfbtemp``

   :Units: Kelvin
   :Description: ice temperature at the basal surface of ice shelves
   :Standard name: ---

#. ``ssa_bc_mask``

   :Units: ---
   :Description: Dirichlet boundary mask
   :Standard name: ---

#. ``ssa_bc_vel``

   - ``u_ssa_bc``

     :Units: m year-1
     :Description: X-component of the SSA velocity boundary conditions
     :Standard name: ---

   - ``v_ssa_bc``

     :Units: m year-1
     :Description: Y-component of the SSA velocity boundary conditions
     :Standard name: ---

#. ``strain_rates``

   - ``eigen1``

     :Units: s-1
     :Description: first eigenvalue of the horizontal, vertically-integrated strain rate tensor
     :Standard name: ---

   - ``eigen2``

     :Units: s-1
     :Description: second eigenvalue of the horizontal, vertically-integrated strain rate tensor
     :Standard name: ---

#. ``strainheat``

   :Units: mW m-3
   :Description: rate of strain heating in ice (dissipation heating)
   :Standard name: ---

#. ``surface_layer_mass``

   :Units: kg
   :Description: mass of the surface layer (snow and firn)
   :Standard name: ---

#. ``surface_layer_thickness``

   :Units: meters
   :Description: thickness of the surface layer (snow and firn)
   :Standard name: ---

#. ``taub``

   - ``taub_x``

     :Units: Pa
     :Description: X-component of the shear stress at the base of ice
     :Standard name: ---
     :Comment: this field is purely diagnostic (not used by the model)

   - ``taub_y``

     :Units: Pa
     :Description: Y-component of the shear stress at the base of ice
     :Standard name: ---
     :Comment: this field is purely diagnostic (not used by the model)

#. ``taub_mag``

   :Units: Pa
   :Description: magnitude of the basal shear stress at the base of ice
   :Standard name: ``magnitude_of_land_ice_basal_drag``
   :Comment: this field is purely diagnostic (not used by the model)

#. ``taud``

   - ``taud_x``

     :Units: Pa
     :Description: X-component of the driving shear stress at the base of ice
     :Standard name: ---
     :Comment: this field is purely diagnostic (not used by the model)

   - ``taud_y``

     :Units: Pa
     :Description: Y-component of the driving shear stress at the base of ice
     :Standard name: ---
     :Comment: this field is purely diagnostic (not used by the model)

#. ``taud_mag``

   :Units: Pa
   :Description: magnitude of the gravitational driving stress at the base of ice
   :Standard name: ---
   :Comment: this field is purely diagnostic (not used by the model)

#. ``tauxz``

   :Units: Pa
   :Description: shear stress xz component (in shallow ice approximation SIA)
   :Standard name: ---

#. ``tauyz``

   :Units: Pa
   :Description: shear stress yz component (in shallow ice approximation SIA)
   :Standard name: ---

#. ``temp``

   :Units: K
   :Description: ice temperature
   :Standard name: ``land_ice_temperature``

#. ``temp_pa``

   :Units: deg_C
   :Description: pressure-adjusted ice temperature (degrees above pressure-melting point)
   :Standard name: ---

#. ``tempbase``

   :Units: K
   :Description: ice temperature at the base of ice
   :Standard name: ``land_ice_basal_temperature``

#. ``tempicethk``

   :Units: m
   :Description: temperate ice thickness (total column content)
   :Standard name: ---

#. ``tempicethk_basal``

   :Units: m
   :Description: thickness of the basal layer of temperate ice
   :Standard name: ---

#. ``temppabase``

   :Units: Celsius
   :Description: pressure-adjusted ice temperature at the base of ice
   :Standard name: ---

#. ``tempsurf``

   :Units: K
   :Description: ice temperature at 1m below the ice surface
   :Standard name: ``temperature_at_ground_level_in_snow_or_firn``

#. ``tendency_of_ice_amount``

   :Units: kg m-2 year-1
   :Description: rate of change of the ice amount
   :Standard name: ---

#. ``tendency_of_ice_amount_due_to_basal_mass_flux``

   :Units: kg m-2 year-1
   :Description: average basal mass flux over reporting interval
   :Standard name: ---
   :Comment: positive flux corresponds to ice gain

#. ``tendency_of_ice_amount_due_to_conservation_error``

   :Units: kg m-2 year-1
   :Description: average mass conservation error flux over reporting interval
   :Standard name: ---
   :Comment: positive flux corresponds to ice gain

#. ``tendency_of_ice_amount_due_to_discharge``

   :Units: kg m-2 year-1
   :Description: discharge (calving and frontal melt) flux
   :Standard name: ``land_ice_specific_mass_flux_due_to_calving_and_ice_front_melting``
   :Comment: positive flux corresponds to ice gain

#. ``tendency_of_ice_amount_due_to_flow``

   :Units: kg m-2 year-1
   :Description: rate of change of ice amount due to flow
   :Standard name: ---
   :Comment: positive flux corresponds to ice gain

#. ``tendency_of_ice_amount_due_to_surface_mass_flux``

   :Units: kg m-2 year-1
   :Description: average surface mass flux over reporting interval
   :Standard name: ---
   :Comment: positive flux corresponds to ice gain

#. ``tendency_of_ice_mass``

   :Units: Gt year-1
   :Description: rate of change of the ice mass
   :Standard name: ---

#. ``tendency_of_ice_mass_due_to_basal_mass_flux``

   :Units: Gt year-1
   :Description: average basal mass flux over reporting interval
   :Standard name: ``tendency_of_land_ice_mass_due_to_basal_mass_balance``
   :Comment: positive flux corresponds to ice gain

#. ``tendency_of_ice_mass_due_to_conservation_error``

   :Units: Gt year-1
   :Description: average mass conservation error flux over reporting interval
   :Standard name: ---
   :Comment: positive flux corresponds to ice gain

#. ``tendency_of_ice_mass_due_to_discharge``

   :Units: Gt year-1
   :Description: discharge (calving and frontal melt) flux
   :Standard name: ---
   :Comment: positive flux corresponds to ice gain

#. ``tendency_of_ice_mass_due_to_flow``

   :Units: Gt year-1
   :Description: rate of change of ice mass due to flow
   :Standard name: ---
   :Comment: positive flux corresponds to ice gain

#. ``tendency_of_ice_mass_due_to_surface_mass_flux``

   :Units: Gt year-1
   :Description: average surface mass flux over reporting interval
   :Standard name: ---
   :Comment: positive flux corresponds to ice gain

#. ``thk``

   :Units: m
   :Description: land ice thickness
   :Standard name: ``land_ice_thickness``

#. ``thksmooth``

   :Units: m
   :Description: thickness relative to smoothed bed elevation in Schoof's (2003) theory of bed roughness in SIA
   :Standard name: ---

#. ``tillwat``

   :Units: m
   :Description: effective thickness of subglacial water stored in till
   :Standard name: ---

#. ``topg``

   :Units: m
   :Description: bedrock surface elevation
   :Standard name: ``bedrock_altitude``

#. ``topg_sl_adjusted``

   :Units: meters
   :Description: sea-level adjusted bed topography (zero is at sea level)
   :Standard name: ---

#. ``topgsmooth``

   :Units: m
   :Description: smoothed bed elevation in Schoof's (2003) theory of bed roughness in SIA
   :Standard name: ---

#. ``usurf``

   :Units: m
   :Description: ice upper surface elevation
   :Standard name: ``surface_altitude``

#. ``uvel``

   :Units: m year-1
   :Description: horizontal velocity of ice in the X direction
   :Standard name: ``land_ice_x_velocity``

#. ``velbar``

   - ``ubar``

     :Units: m year-1
     :Description: vertical mean of horizontal ice velocity in the X direction
     :Standard name: ``land_ice_vertical_mean_x_velocity``

   - ``vbar``

     :Units: m year-1
     :Description: vertical mean of horizontal ice velocity in the Y direction
     :Standard name: ``land_ice_vertical_mean_y_velocity``

#. ``velbar_mag``

   :Units: m year-1
   :Description: magnitude of vertically-integrated horizontal velocity of ice
   :Standard name: ---

#. ``velbase``

   - ``uvelbase``

     :Units: m year-1
     :Description: x-component of the horizontal velocity of ice at the base of ice
     :Standard name: ``land_ice_basal_x_velocity``

   - ``vvelbase``

     :Units: m year-1
     :Description: y-component of the horizontal velocity of ice at the base of ice
     :Standard name: ``land_ice_basal_y_velocity``

#. ``velbase_mag``

   :Units: m year-1
   :Description: magnitude of horizontal velocity of ice at base of ice
   :Standard name: ---

#. ``velsurf``

   - ``uvelsurf``

     :Units: m year-1
     :Description: x-component of the horizontal velocity of ice at ice surface
     :Standard name: ``land_ice_surface_x_velocity``

   - ``vvelsurf``

     :Units: m year-1
     :Description: y-component of the horizontal velocity of ice at ice surface
     :Standard name: ``land_ice_surface_y_velocity``

#. ``velsurf_mag``

   :Units: m year-1
   :Description: magnitude of horizontal velocity of ice at ice surface
   :Standard name: ---

#. ``vonmises_stress``

   :Units: Pascal
   :Description: tensile von Mises stress
   :Standard name: ---

#. ``vvel``

   :Units: m year-1
   :Description: horizontal velocity of ice in the Y direction
   :Standard name: ``land_ice_y_velocity``

#. ``wallmelt``

   :Units: m year-1
   :Description: wall melt into subglacial hydrology layer from (turbulent) dissipation of energy in transportable water
   :Standard name: ---

#. ``wvel``

   :Units: m year-1
   :Description: vertical velocity of ice, relative to geoid
   :Standard name: ---

#. ``wvel_rel``

   :Units: m year-1
   :Description: vertical velocity of ice, relative to base of ice directly below
   :Standard name: ---

#. ``wvelbase``

   :Units: m year-1
   :Description: vertical velocity of ice at the base of ice, relative to the geoid
   :Standard name: ``land_ice_basal_upward_velocity``

#. ``wvelsurf``

   :Units: m year-1
   :Description: vertical velocity of ice at ice surface, relative to the geoid
   :Standard name: ``land_ice_surface_upward_velocity``

.. _sec-ts_vars:

Scalar time-series
------------------

#. ``area_glacierized``

   :Units: m2
   :Description: glacierized area
   :Standard name: ---

#. ``area_glacierized_cold_base``

   :Units: m2
   :Description: glacierized area where basal ice is cold
   :Standard name: ---

#. ``area_glacierized_floating``

   :Units: m2
   :Description: area of ice shelves in glacierized areas
   :Standard name: ---

#. ``area_glacierized_grounded``

   :Units: m2
   :Description: area of grounded ice in glacierized areas
   :Standard name: ---

#. ``area_glacierized_temperate_base``

   :Units: m2
   :Description: glacierized area where basal ice is temperate
   :Standard name: ---

#. ``basal_mass_flux_floating``

   :Units: kg year-1
   :Description: total sub-shelf ice flux
   :Standard name: ---
   :Comment: positive means ice gain

#. ``basal_mass_flux_grounded``

   :Units: kg year-1
   :Description: total over grounded ice domain of basal mass flux
   :Standard name: ---
   :Comment: positive means ice gain

#. ``dt``

   :Units: year
   :Description: mass continuity time step
   :Standard name: ---

#. ``enthalpy_glacierized``

   :Units: J
   :Description: enthalpy of the ice in glacierized areas
   :Standard name: ---

#. ``enthalpy_nonglacierized``

   :Units: J
   :Description: enthalpy of the ice, including seasonal cover
   :Standard name: ---

#. ``limnsw``

   :Units: kg
   :Description: mass of the ice not displacing sea water
   :Standard name: ---

#. ``liquified_ice_flux``

   :Units: m3 / year
   :Description: rate of ice loss due to liquefaction, averaged over the reporting interval
   :Standard name: ---
   :Comment: positive means ice loss

#. ``mass_glacierized``

   :Units: kg
   :Description: mass of the ice in glacierized areas
   :Standard name: ---

#. ``mass_nonglacierized``

   :Units: kg
   :Description: mass of the ice, including seasonal cover
   :Standard name: ---

#. ``mass_rate_of_change_glacierized``

   :Units: kg year-1
   :Description: rate of change of the mass of ice in glacierized areas
   :Standard name: ---

#. ``mass_rate_of_change_nonglacierized``

   :Units: kg year-1
   :Description: rate of change of the mass of ice, including seasonal cover
   :Standard name: ---

#. ``max_diffusivity``

   :Units: m2 s-1
   :Description: maximum diffusivity
   :Standard name: ---

#. ``max_hor_vel``

   :Units: m year-1
   :Description: maximum abs component of horizontal ice velocity over grid in last time step during time-series reporting interval
   :Standard name: ---

#. ``slvol``

   :Units: m
   :Description: total sea-level relevant ice IN SEA-LEVEL EQUIVALENT
   :Standard name: ---

#. ``tendency_of_ice_mass``

   :Units: kg year-1
   :Description: rate of change of the mass of ice, including seasonal cover
   :Standard name: ---

#. ``tendency_of_ice_mass_due_to_basal_mass_balance``

   :Units: kg year-1
   :Description: total over ice domain of bottom surface ice mass flux
   :Standard name: ---
   :Comment: positive means ice gain

#. ``tendency_of_ice_mass_due_to_conservation_error``

   :Units: kg year-1
   :Description: total numerical flux needed to preserve non-negativity of ice thickness
   :Standard name: ---
   :Comment: positive means ice gain

#. ``tendency_of_ice_mass_due_to_discharge``

   :Units: kg year-1
   :Description: discharge (calving & icebergs) flux
   :Standard name: ---
   :Comment: positive means ice gain

#. ``tendency_of_ice_mass_due_to_influx``

   :Units: kg year-1
   :Description: rate of change of the mass of ice due to influx (i.e. prescribed ice thickness)
   :Standard name: ---

#. ``tendency_of_ice_mass_due_to_surface_mass_balance``

   :Units: kg year-1
   :Description: total over ice domain of top surface ice mass flux
   :Standard name: ---
   :Comment: positive means ice gain

#. ``volume_glacierized``

   :Units: m3
   :Description: volume of the ice in glacierized areas
   :Standard name: ---

#. ``volume_glacierized_cold``

   :Units: m3
   :Description: volume of cold ice in glacierized areas
   :Standard name: ---

#. ``volume_glacierized_floating``

   :Units: m3
   :Description: volume of ice shelves in glacierized areas
   :Standard name: ---

#. ``volume_glacierized_grounded``

   :Units: m3
   :Description: volume of grounded ice in glacierized areas
   :Standard name: ---

#. ``volume_glacierized_temperate``

   :Units: m3
   :Description: volume of temperate ice in glacierized areas
   :Standard name: ---

#. ``volume_nonglacierized``

   :Units: m3
   :Description: volume of the ice, including seasonal cover
   :Standard name: ---

#. ``volume_nonglacierized_cold``

   :Units: m3
   :Description: volume of cold ice, including seasonal cover
   :Standard name: ---

#. ``volume_nonglacierized_temperate``

   :Units: m3
   :Description: volume of temperate ice, including seasonal cover
   :Standard name: ---

#. ``volume_rate_of_change_glacierized``

   :Units: m3 year-1
   :Description: rate of change of the ice volume in glacierized areas
   :Standard name: ---

#. ``volume_rate_of_change_nonglacierized``

   :Units: m3 year-1
   :Description: rate of change of the ice volume, including seasonal cover
   :Standard name: ---
