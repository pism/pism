
.. DO NOT EDIT: This file was automatically generated using config_parameters.py. Edit src/pism_config.cdl instead.

.. include:: ../../global.rst

.. _sec-parameter-list:

List of configuration parameters
================================

Each parameter can be set using the command-line option consisting of a dash followed by
the parameter name. For example,

.. code-block:: none

   -constants.standard_gravity 10

sets the acceleration due to gravity (parameter :config:`constants.standard_gravity`) to
`10`.



#. :config:`age.enabled`

   :Value: no
   :Option: :opt:`-age`
   :Description: Solve age equation (advection equation for ice age).

#. :config:`age.initial_value`

   :Value: 0 (years)
   :Option: *no short option*
   :Description: Initial age of ice

#. :config:`atmosphere.fausto_air_temp.c_ma`

   :Value: -0.718900 (Kelvin / degree_north)
   :Option: *no short option*
   :Description: latitude-dependence coefficient for formula (1) in :cite:`Faustoetal2009`

#. :config:`atmosphere.fausto_air_temp.c_mj`

   :Value: -0.158500 (Kelvin / degree_north)
   :Option: *no short option*
   :Description: latitude-dependence coefficient for formula (2) in :cite:`Faustoetal2009`

#. :config:`atmosphere.fausto_air_temp.d_ma`

   :Value: 314.980000 (Kelvin)
   :Option: *no short option*
   :Description: 41.83+273.15; base temperature for formula (1) in :cite:`Faustoetal2009`

#. :config:`atmosphere.fausto_air_temp.d_mj`

   :Value: 287.850000 (Kelvin)
   :Option: *no short option*
   :Description: = 14.70+273.15; base temperature for formula (2) in :cite:`Faustoetal2009`

#. :config:`atmosphere.fausto_air_temp.gamma_ma`

   :Value: -0.006309 (Kelvin / meter)
   :Option: *no short option*
   :Description: = -6.309 / 1km; mean slope lapse rate for formula (1) in :cite:`Faustoetal2009`

#. :config:`atmosphere.fausto_air_temp.gamma_mj`

   :Value: -0.005426 (Kelvin / meter)
   :Option: *no short option*
   :Description: = -5.426 / 1km; mean slope lapse rate for formula (2) in :cite:`Faustoetal2009`

#. :config:`atmosphere.fausto_air_temp.kappa_ma`

   :Value: 0.067200 (Kelvin / degree_west)
   :Option: *no short option*
   :Description: longitude-dependence coefficient for formula (1) in :cite:`Faustoetal2009`

#. :config:`atmosphere.fausto_air_temp.kappa_mj`

   :Value: 0.051800 (Kelvin / degree_west)
   :Option: *no short option*
   :Description: longitude-dependence coefficient for formula (2) in :cite:`Faustoetal2009`

#. :config:`atmosphere.fausto_air_temp.summer_peak_day`

   :Value: 196
   :Option: *no short option*
   :Description: day of year for July 15; used in corrected formula (4) in :cite:`Faustoetal2009`

#. :config:`atmosphere.precip_exponential_factor_for_temperature`

   :Value: 0.070417 (Kelvin-1)
   :Option: *no short option*
   :Description: = 0.169/2.4; in SeaRISE-Greenland formula for paleo-precipitation from present; a 7.3% change of precipitation rate for every one degC of temperature change :cite:`Huybrechts02`

#. :config:`basal_resistance.beta_ice_free_bedrock`

   :Value: 1.800000e+09 (Pascal second meter-1)
   :Option: *no short option*
   :Description: value is for ice stream E from :cite:`HulbeMacAyeal`; thus sliding velocity, but we hope it doesn't matter much; at 100 m/year the linear sliding law gives 57040 Pa basal shear stress

#. :config:`basal_resistance.plastic.regularization`

   :Value: 0.010000 (meter / year)
   :Option: :opt:`-plastic_reg`
   :Description: Set the value of `\epsilon` regularization of plastic till; this is the second `\epsilon` in formula (4.1) in :cite:`SchoofStream`

#. :config:`basal_resistance.pseudo_plastic.enabled`

   :Value: no
   :Option: :opt:`-pseudo_plastic`
   :Description: Use the pseudo-plastic till model (basal sliding law).

#. :config:`basal_resistance.pseudo_plastic.q`

   :Value: 0.250000 (pure number)
   :Option: :opt:`-pseudo_plastic_q`
   :Description: The exponent of the pseudo-plastic basal resistance model

#. :config:`basal_resistance.pseudo_plastic.sliding_scale_factor`

   :Value: -1 (1)
   :Option: :opt:`-sliding_scale_factor_reduces_tauc`
   :Description: divides pseudo-plastic tauc (yield stress) by given factor; this would increase sliding by given factor in absence of membrane stresses; not used if negative or zero; not used by default

#. :config:`basal_resistance.pseudo_plastic.u_threshold`

   :Value: 100 (meter / year)
   :Option: :opt:`-pseudo_plastic_uthreshold`
   :Description: threshold velocity of the pseudo-plastic sliding law

#. :config:`basal_yield_stress.add_transportable_water`

   :Value: no
   :Option: :opt:`-tauc_add_transportable_water`
   :Description: If 'yes' then the water amount in the transport system is added to tillwat in determining tauc (in the Mohr-Coulomb relation).  Normally only the water in the till is used.

#. :config:`basal_yield_stress.constant.value`

   :Value: 200000 (Pascal)
   :Option: :opt:`-tauc`
   :Description: fill value for yield stress for basal till (plastic or pseudo-plastic model); note 2 x 10^5 Pa = 2.0 bar is quite strong and little sliding should occur without an explicit tauc choice altering this default

#. :config:`basal_yield_stress.ice_free_bedrock`

   :Value: 1000000 (Pascal)
   :Option: :opt:`-high_tauc`
   :Description: the 'high' yield stress value used in grounded ice-free areas.

#. :config:`basal_yield_stress.model`

   :Value: mohr_coulomb
   :Choices: ``constant, mohr_coulomb``
   :Option: :opt:`-yield_stress`
   :Description: The basal yield stress model to use when sliding is active.

#. :config:`basal_yield_stress.mohr_coulomb.till_cohesion`

   :Value: 0 (Pascal)
   :Option: :opt:`-till_cohesion`
   :Description: cohesion of till; = c_0 in most references; note Schoof uses zero but Paterson pp 168--169 gives range 0--40 kPa; but Paterson notes that '... all the pairs c_0 and phi in the table would give a yield stress for Ice Stream B that exceeds the basal shear stress there...'

#. :config:`basal_yield_stress.mohr_coulomb.till_compressibility_coefficient`

   :Value: 0.120000 (pure number)
   :Option: :opt:`-till_compressibility_coefficient`
   :Description: coefficient of compressiblity of till; value from :cite:`Tulaczyketal2000`

#. :config:`basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden`

   :Value: 0.020000 (pure number)
   :Option: :opt:`-till_effective_fraction_overburden`
   :Description: `\delta` in notes; `N_0 = \delta P_o` where `P_o` is overburden pressure; `N_0` is reference (low) value of effective pressure (i.e. normal stress); `N_0` scales with overburden pressure unlike :cite:`Tulaczyketal2000`; default value from Greenland and Antarctic model runs

#. :config:`basal_yield_stress.mohr_coulomb.till_log_factor_transportable_water`

   :Value: 0.100000 (meters)
   :Option: :opt:`-till_log_factor_transportable_water`
   :Description: If basal_yield_stress.add_transportable_water = yes then the water amount in the transport system is added to tillwat in determining tauc.  Normally only the water in the till is used.  This factor multiplies the logarithm in that formula.

#. :config:`basal_yield_stress.mohr_coulomb.till_phi_default`

   :Value: 30 (degrees)
   :Option: :opt:`-plastic_phi`
   :Description: fill value for till friction angle

#. :config:`basal_yield_stress.mohr_coulomb.till_reference_effective_pressure`

   :Value: 1000 (Pascal)
   :Option: *no short option*
   :Description: reference effective pressure N_0; value from :cite:`Tulaczyketal2000`

#. :config:`basal_yield_stress.mohr_coulomb.till_reference_void_ratio`

   :Value: 0.690000 (pure number)
   :Option: :opt:`-till_reference_void_ratio`
   :Description: void ratio at reference effective pressure N_0; value from :cite:`Tulaczyketal2000`

#. :config:`basal_yield_stress.mohr_coulomb.topg_to_phi.enabled`

   :Value: no
   :Option: *no short option*
   :Description: If THE OPTION -topg_to_phi IS SET THEN THIS WILL BE SET TO 'yes', and then MohrCoulombYieldStress will initialize the tillphi field using a piece-wise linear function of depth described by four parameters.

#. :config:`basal_yield_stress.mohr_coulomb.topg_to_phi.phi_max`

   :Value: 15 (degrees)
   :Option: *no short option*
   :Description: upper value of the till friction angle; see the implementation of MohrCoulombYieldStress

#. :config:`basal_yield_stress.mohr_coulomb.topg_to_phi.phi_min`

   :Value: 5 (degrees)
   :Option: *no short option*
   :Description: lower value of the till friction angle; see the implementation of MohrCoulombYieldStress

#. :config:`basal_yield_stress.mohr_coulomb.topg_to_phi.topg_max`

   :Value: 1000 (meters)
   :Option: *no short option*
   :Description: the elevation at which the upper value of the till friction angle is used; see the implementation of MohrCoulombYieldStress

#. :config:`basal_yield_stress.mohr_coulomb.topg_to_phi.topg_min`

   :Value: -1000 (meters)
   :Option: *no short option*
   :Description: the elevation at which the lower value of the till friction angle is used; see the implementation of MohrCoulombYieldStress

#. :config:`basal_yield_stress.slippery_grounding_lines`

   :Value: no
   :Option: :opt:`-tauc_slippery_grounding_lines`
   :Description: If yes, at icy grounded locations with bed elevations below sea level, within one cell of floating ice or ice-free ocean, make tauc as low as possible from the Mohr-Coulomb relation.  Specifically, at such locations replace the normally-computed tauc from the Mohr-Coulomb relation, which uses the effective pressure from the modeled amount of water in the till, by the minimum value of tauc from Mohr-Coulomb, i.e. by using the effective pressure corresponding to the maximum amount of till-stored water.  Does not alter the modeled or reported amount of till water, nor does this mechanism affect water conservation.

#. :config:`bed_deformation.lc.elastic_model`

   :Value: no
   :Option: :opt:`-bed_def_lc_elastic_model`
   :Description: Use the elastic part of the Lingle-Clark bed deformation model.

#. :config:`bed_deformation.mantle_density`

   :Value: 3300 (kg meter-3)
   :Option: *no short option*
   :Description: half-space (mantle) density used by the bed deformation model. See :cite:`LingleClark`, :cite:`BLKfastearth`

#. :config:`bed_deformation.lithosphere_flexural_rigidity`

   :Value: 5.000000e+24 (Newton meter)
   :Option: *no short option*
   :Description: lithosphere flexural rigidity used by the bed deformation model. See :cite:`LingleClark`, :cite:`BLKfastearth`

#. :config:`bed_deformation.mantle_viscosity`

   :Value: 1.000000e+21 (Pascal second)
   :Option: *no short option*
   :Description: half-space (mantle) viscosity used by the bed deformation model. See :cite:`LingleClark`, :cite:`BLKfastearth`

#. :config:`bed_deformation.model`

   :Value: none
   :Choices: ``none, iso, lc``
   :Option: :opt:`-bed_def`
   :Description: Selects a bed deformation model to use. 'iso' is point-wise isostasy, 'lc' is the Lingle-Clark model (see :cite:`LingleClark`, requires FFTW3).

#. :config:`bed_deformation.update_interval`

   :Value: 10 (years)
   :Option: *no short option*
   :Description: Interval between bed deformation updates

#. :config:`bed_deformation.bed_topography_delta_file`

   :Value: *no default*
   :Option: :opt:`-topg_delta_file`
   :Description: The name of the file to read the topg_delta from. This field is added to the bed topography during initialization.

#. :config:`bed_deformation.bed_uplift_file`

   :Value: *no default*
   :Option: :opt:`-uplift_file`
   :Description: The name of the file to read the uplift (dbdt) from. Leave empty to read it from an input file or a regridding file.

#. :config:`bootstrapping.defaults.bed`

   :Value: 1 (meters)
   :Option: *no short option*
   :Description: bed elevation value to use if topg (bedrock_altitude) variable is absent in bootstrapping file

#. :config:`bootstrapping.defaults.bmelt`

   :Value: 0 (meter / second)
   :Option: *no short option*
   :Description: basal melt rate value to use if variable bmelt is absent in bootstrapping file

#. :config:`bootstrapping.defaults.bwat`

   :Value: 0 (meters)
   :Option: *no short option*
   :Description: till water thickness value to use if variable tillwat is absent in bootstrapping file

#. :config:`bootstrapping.defaults.bwp`

   :Value: 0 (Pascal)
   :Option: *no short option*
   :Description: basal water pressure value to use if variable bwp is absent in bootstrapping file; most hydrology models do not use this value because bwp is diagnostic

#. :config:`bootstrapping.defaults.enwat`

   :Value: 0 (meters)
   :Option: *no short option*
   :Description: effective englacial water thickness value to use if variable enwat is absent in bootstrapping file

#. :config:`bootstrapping.defaults.geothermal_flux`

   :Value: 0.042000 (Watt meter-2)
   :Option: *no short option*
   :Description: geothermal flux value to use if bheatflx variable is absent in bootstrapping file

#. :config:`bootstrapping.defaults.ice_thickness`

   :Value: 0 (meters)
   :Option: *no short option*
   :Description: thickness value to use if thk (land_ice_thickness) variable is absent in bootstrapping file

#. :config:`bootstrapping.defaults.tillwat`

   :Value: 0 (meters)
   :Option: *no short option*
   :Description: till water thickness value to use if variable tillwat is absent in bootstrapping file

#. :config:`bootstrapping.defaults.uplift`

   :Value: 0 (meter / second)
   :Option: *no short option*
   :Description: uplift value to use if dbdt variable is absent in bootstrapping file

#. :config:`bootstrapping.temperature_heuristic`

   :Value: smb
   :Choices: ``smb, quartic_guess``
   :Option: :opt:`-boot_temperature_heuristic`
   :Description: The heuristic to use to initialize ice temperature during bootstrapping: 'sbm' uses the new method using the surface mass balance, surface temperature, and the geothermal flux, 'quartic_guess' uses the old method using the surface temperature and the geothermal flux.

#. :config:`calving.eigen_calving.K`

   :Value: 0 (meter second)
   :Option: :opt:`-eigen_calving_K`
   :Description: Set proportionality constant to determine calving rate from strain rates.  Note references :cite:`Levermannetal2012`, :cite:`Martinetal2011` use K in range 10^9 to 3 x 10^11 m a, that is, 3 x 10^16 to 10^19 m s.

#. :config:`calving.float_kill.margin_only`

   :Value: no
   :Option: :opt:`-float_kill_margin_only`
   :Description: Apply float_kill at ice margin cells only.

#. :config:`calving.float_kill.calve_near_grounding_line`

   :Value: yes
   :Option: :opt:`-float_kill_calve_near_grounding_line`
   :Description: Calve floating ice near the grounding line.

#. :config:`calving.front_retreat.wrap_around`

   :Value: false
   :Option: :opt:`-calving_wrap_around`
   :Description: If true, wrap around domain boundaries. This may be needed in some regional synthetic geometry setups.

#. :config:`calving.front_retreat.use_cfl`

   :Value: false
   :Option: :opt:`-calving_cfl`
   :Description: apply CFL criterion for eigen-calving rate front retreat

#. :config:`calving.methods`

   :Value: *no default*
   :Option: :opt:`-calving`
   :Description: comma-separated list of calving methods; one or more of 'eigen_calving', 'ocean_kill', 'float_kill', 'thickness_calving'

#. :config:`calving.thickness_calving.threshold`

   :Value: 50 (meters)
   :Option: :opt:`-thickness_calving_threshold`
   :Description: When terminal ice thickness of floating ice shelf is less than this threshold, it will be calved off.

#. :config:`calving.thickness_calving.threshold_file`

   :Value: *no default*
   :Option: :opt:`-thickness_calving_threshold_file`
   :Description: Name of the file containing the spatially-variable thickness calving threshold. 

#. :config:`calving.vonmises.sigma_max`

   :Value: 1000000 (Pa)
   :Option: :opt:`-vonmises_calving_sigma_max`
   :Description: Set maximum tensile stress.  Note references :cite:`Morlighem2016` use 1.0e6 Pa.

#. :config:`climate_forcing.buffer_size`

   :Value: 60
   :Option: *no short option*
   :Description: number of 2D climate forcing records to keep in memory; = 5 years of monthly records

#. :config:`climate_forcing.evaluations_per_year`

   :Value: 52
   :Option: *no short option*
   :Description: length of the time-series used to compute temporal averages of forcing data (such as mean annual temperature)

#. :config:`constants.fresh_water.density`

   :Value: 1000 (kg meter-3)
   :Option: *no short option*
   :Description: density of fresh water

#. :config:`constants.fresh_water.latent_heat_of_fusion`

   :Value: 334000 (Joule / kg)
   :Option: *no short option*
   :Description: latent heat of fusion for water :cite:`AschwandenBlatter`

#. :config:`constants.fresh_water.melting_point_temperature`

   :Value: 273.150000 (Kelvin)
   :Option: *no short option*
   :Description: melting point of pure water

#. :config:`constants.fresh_water.specific_heat_capacity`

   :Value: 4170 (Joule / (kg Kelvin))
   :Option: *no short option*
   :Description: at melting point T_0 :cite:`AschwandenBlatter`

#. :config:`constants.ice.beta_Clausius_Clapeyron`

   :Value: 7.900000e-08 (Kelvin / Pascal)
   :Option: *no short option*
   :Description: Clausius-Clapeyron constant relating melting temperature and pressure: `\beta = dT / dP` :cite:`Luethi2002`

#. :config:`constants.ice.density`

   :Value: 910 (kg meter-3)
   :Option: *no short option*
   :Description: `\rho_i`; density of ice in ice sheet

#. :config:`constants.ice.grain_size`

   :Value: 1 (mm)
   :Option: :opt:`-ice_grain_size`
   :Description: Default constant ice grain size to use with the Goldsby-Kohlstedt :cite:`GoldsbyKohlstedt` flow law

#. :config:`constants.ice.specific_heat_capacity`

   :Value: 2009 (Joule / (kg Kelvin))
   :Option: *no short option*
   :Description: specific heat capacity of pure ice at melting point T_0

#. :config:`constants.ice.thermal_conductivity`

   :Value: 2.100000 (Joule / (meter Kelvin second))
   :Option: *no short option*
   :Description: = W m-1 K-1; thermal conductivity of pure ice

#. :config:`constants.ideal_gas_constant`

   :Value: 8.314410 (Joule / (mol Kelvin))
   :Option: *no short option*
   :Description: ideal gas constant

#. :config:`constants.sea_water.density`

   :Value: 1028 (kg meter-3)
   :Option: *no short option*
   :Description: density of sea water

#. :config:`constants.sea_water.specific_heat_capacity`

   :Value: 3985 (Joule / (kg Kelvin))
   :Option: *no short option*
   :Description: at 35 psu, value taken from http://www.kayelaby.npl.co.uk/general_physics/2_7/2_7_9.html

#. :config:`constants.standard_gravity`

   :Value: 9.810000 (meter second-2)
   :Option: *no short option*
   :Description: acceleration due to gravity on Earth geoid

#. :config:`energy.allow_temperature_above_melting`

   :Value: no
   :Option: *no short option*
   :Description: If set to 'yes', allow temperatures above the pressure-malting point in the cold mode temperature code. Used by some verifiaction tests.

#. :config:`energy.basal_melt.use_grounded_cell_fraction`

   :Value: true
   :Option: :opt:`-subgl_basal_melt`
   :Description: If geometry.grounded_cell_fraction is set, use the fractional floatation mask to interpolate the basal melt rate at the grounding line between grounded and floating values.

#. :config:`energy.bedrock_thermal_conductivity`

   :Value: 3 (Joule / (meter Kelvin second))
   :Option: *no short option*
   :Description: = W m-1 K-1; for bedrock used in thermal model :cite:`RitzEISMINT`

#. :config:`energy.bedrock_thermal_density`

   :Value: 3300 (kg meter-3)
   :Option: *no short option*
   :Description: for bedrock used in thermal model

#. :config:`energy.bedrock_thermal_specific_heat_capacity`

   :Value: 1000 (Joule / (kg Kelvin))
   :Option: *no short option*
   :Description: for bedrock used in thermal model :cite:`RitzEISMINT`

#. :config:`energy.drainage_maximum_rate`

   :Value: 1.584438e-09 (second-1)
   :Option: *no short option*
   :Description: 0.05 year-1; maximum rate at which liquid water fraction in temperate ice could possibly drain; see :cite:`AschwandenBuelerKhroulevBlatter`

#. :config:`energy.drainage_target_water_fraction`

   :Value: 0.010000 (1)
   :Option: *no short option*
   :Description: liquid water fraction (omega) above which drainage occurs, but below which there is no drainage; see :cite:`AschwandenBuelerKhroulevBlatter`

#. :config:`energy.enabled`

   :Value: yes
   :Option: *no short option*
   :Description: Solve energy conservation equations.

#. :config:`energy.enthalpy_cold_bulge_max`

   :Value: 60270 (Joule / kg)
   :Option: *no short option*
   :Description: = (2009 J kg-1 K-1) * (30 K); maximum amount by which advection can reduce the enthalpy of a column of ice below its surface enthalpy value

#. :config:`energy.max_low_temperature_count`

   :Value: 10
   :Option: :opt:`-max_low_temps`
   :Description: Maximum number of grid points with ice temperature below energy.minimum_allowed_temperature.

#. :config:`energy.minimum_allowed_temperature`

   :Value: 200 (Kelvin)
   :Option: :opt:`-low_temp`
   :Description: Minimum allowed ice temperature

#. :config:`energy.temperate_ice_enthalpy_conductivity_ratio`

   :Value: 0.100000 (pure number)
   :Option: *no short option*
   :Description: K in cold ice is multiplied by this fraction to give K0 in :cite:`AschwandenBuelerKhroulevBlatter`

#. :config:`energy.temperature_based`

   :Value: no
   :Option: *no short option*
   :Description: Use cold ice (i.e. not polythermal) methods.

#. :config:`energy.temperature_dependent_thermal_conductivity`

   :Value: no
   :Option: :opt:`-vark`
   :Description: If yes, use varkenthSystemCtx class in the energy step. It is base on formula (4.37) in :cite:`GreveBlatter2009`. Otherwise use enthSystemCtx, which has temperature-independent thermal conductivity set by constant ice.thermal_conductivity.

#. :config:`enthalpy_converter.T_reference`

   :Value: 223.150000 (Kelvin)
   :Option: *no short option*
   :Description: = T_0 in enthalpy formulas in :cite:`AschwandenBuelerKhroulevBlatter`

#. :config:`enthalpy_converter.relaxed_is_temperate_tolerance`

   :Value: 0.001000 (Kelvin)
   :Option: *no short option*
   :Description: Tolerance within which ice is treated as temperate (cold-ice mode and diagnostics).

#. :config:`flow_law.Hooke.A`

   :Value: 4.421650e-09 (Pascal-3 second-1)
   :Option: *no short option*
   :Description: `A_{\text{Hooke}} = (1/B_0)^n` where n=3 and B_0 = 1.928 `a^{1/3}` Pa. See :cite:`Hooke`

#. :config:`flow_law.Hooke.C`

   :Value: 0.166120 (Kelvin^{flow_law.Hooke.k})
   :Option: *no short option*
   :Description: See :cite:`Hooke`

#. :config:`flow_law.Hooke.Q`

   :Value: 78800 (Joule / mol)
   :Option: *no short option*
   :Description: Activation energy, see :cite:`Hooke`

#. :config:`flow_law.Hooke.Tr`

   :Value: 273.390000 (Kelvin)
   :Option: *no short option*
   :Description: See :cite:`Hooke`

#. :config:`flow_law.Hooke.k`

   :Value: 1.170000 (pure number)
   :Option: *no short option*
   :Description: See :cite:`Hooke`

#. :config:`flow_law.Paterson_Budd.A_cold`

   :Value: 3.610000e-13 (Pascal-3 / second)
   :Option: *no short option*
   :Description: Paterson-Budd A_cold, see :cite:`PatersonBudd`

#. :config:`flow_law.Paterson_Budd.A_warm`

   :Value: 1730 (Pascal-3 / second)
   :Option: *no short option*
   :Description: Paterson-Budd A_warm, see :cite:`PatersonBudd`

#. :config:`flow_law.Paterson_Budd.Q_cold`

   :Value: 60000 (Joule / mol)
   :Option: *no short option*
   :Description: Paterson-Budd Q_cold, see :cite:`PatersonBudd`

#. :config:`flow_law.Paterson_Budd.Q_warm`

   :Value: 139000 (Joule / mol)
   :Option: *no short option*
   :Description: Paterson-Budd Q_warm, see :cite:`PatersonBudd`

#. :config:`flow_law.Paterson_Budd.T_critical`

   :Value: 263.150000 (Kelvin)
   :Option: *no short option*
   :Description: Paterson-Budd critical temperature, see :cite:`PatersonBudd`

#. :config:`flow_law.Schoof_regularizing_length`

   :Value: 1000 (km)
   :Option: *no short option*
   :Description: Regularizing length (Schoof definition)

#. :config:`flow_law.Schoof_regularizing_velocity`

   :Value: 1 (meter / year)
   :Option: *no short option*
   :Description: Regularizing velocity (Schoof definition)

#. :config:`flow_law.gpbld.water_frac_coeff`

   :Value: 181.250000 (pure number)
   :Option: *no short option*
   :Description: coefficient in Glen-Paterson-Budd flow law for extra dependence of softness on liquid water fraction (omega) :cite:`GreveBlatter2009`, :cite:`LliboutryDuval1985`

#. :config:`flow_law.gpbld.water_frac_observed_limit`

   :Value: 0.010000 (1)
   :Option: *no short option*
   :Description: maximum value of liquid water fraction omega for which softness values are parameterized by :cite:`LliboutryDuval1985`; used in Glen-Paterson-Budd-Lliboutry-Duval flow law; compare :cite:`AschwandenBuelerKhroulevBlatter`

#. :config:`flow_law.isothermal_Glen.ice_softness`

   :Value: 3.168900e-24 (Pascal-3 second-1)
   :Option: *no short option*
   :Description: ice softness used by IsothermalGlenIce :cite:`EISMINT96`

#. :config:`fracture_density.constant_fd`

   :Value: no
   :Option: :opt:`-constant_fd`
   :Description: FIXME

#. :config:`fracture_density.constant_healing`

   :Value: no
   :Option: :opt:`-constant_healing`
   :Description: Constant healing

#. :config:`fracture_density.enabled`

   :Value: no
   :Option: :opt:`-fractures`
   :Description: Calculation of fracture density according to stresses and strain rate field.

#. :config:`fracture_density.fd2d_scheme`

   :Value: no
   :Option: :opt:`-scheme_fd2d`
   :Description: FIXME

#. :config:`fracture_density.fracture_weighted_healing`

   :Value: no
   :Option: :opt:`-fracture_weighted_healing`
   :Description: Fracture weighted healing

#. :config:`fracture_density.include_grounded_ice`

   :Value: no
   :Option: :opt:`-do_frac_on_grounded`
   :Description: model fracture density in grounded areas

#. :config:`fracture_density.lefm`

   :Value: no
   :Option: :opt:`-lefm`
   :Description: FIXME

#. :config:`fracture_density.max_shear_stress`

   :Value: no
   :Option: :opt:`-max_shear`
   :Description: Use the max. shear stress criterion.

#. :config:`fracture_density.phi0`

   :Value: 0 (1)
   :Option: :opt:`-phi0`
   :Description: FIXME

#. :config:`fracture_density.softening_lower_limit`

   :Value: 1 (1)
   :Option: :opt:`-fracture_softening`
   :Description: epsilon in equation (6) in Albrecht and Levermann, 'Fracture-induced softening for large-scale ice dynamics'

#. :config:`fracture_density.write_fields`

   :Value: no
   :Option: :opt:`-write_fd_fields`
   :Description: Writing of fracture density related fields to nc-file.

#. :config:`geometry.grounded_cell_fraction`

   :Value: false
   :Option: :opt:`-subgl`
   :Description: Linear interpolation scheme ('LI' in Gladstone et al. 2010) expanded to two dimensions is used if switched on in order to evaluate the position of the grounding line on a subgrid scale.

#. :config:`geometry.ice_free_thickness_standard`

   :Value: 0.010000 (meters)
   :Option: *no short option*
   :Description: If ice is thinner than this standard then the mask is set to MASK_ICE_FREE_BEDROCK or MASK_ICE_FREE_OCEAN.

#. :config:`geometry.part_grid.enabled`

   :Value: no
   :Option: :opt:`-part_grid`
   :Description: apply partially filled grid cell scheme

#. :config:`geometry.remove_icebergs`

   :Value: no
   :Option: :opt:`-kill_icebergs`
   :Description: identify and kill detached ice-shelf areas

#. :config:`geometry.update.enabled`

   :Value: yes
   :Option: :opt:`-mass`
   :Description: Solve the mass conservation equation

#. :config:`geometry.update.use_basal_melt_rate`

   :Value: yes
   :Option: :opt:`-bmr_in_cont`
   :Description: Include basal melt rate in the continuity equation

#. :config:`grid.allow_extrapolation`

   :Value: no
   :Option: :opt:`-allow_extrapolation`
   :Description: Allow extrapolation during regridding.

#. :config:`grid.Lbz`

   :Value: 0 (meters)
   :Option: *no short option*
   :Description: Thickness of the thermal bedrock layer.

#. :config:`grid.Lx`

   :Value: 1500000 (meters)
   :Option: *no short option*
   :Description: Default computational box is 3000 km x 3000 km (= 2 Lx x 2 Ly) in horizontal.

#. :config:`grid.Ly`

   :Value: 1500000 (meters)
   :Option: *no short option*
   :Description: Default computational box is 3000 km x 3000 km (= 2 Lx x 2 Ly) in horizontal.

#. :config:`grid.Lz`

   :Value: 4000 (meters)
   :Option: *no short option*
   :Description: Height of the computational domain.

#. :config:`grid.Mbz`

   :Value: 1
   :Option: *no short option*
   :Description: Number of thermal bedrock layers; 1 level corresponds to no bedrock.

#. :config:`grid.Mx`

   :Value: 61
   :Option: :opt:`-Mx`
   :Description: Number of grid points in the x direction.

#. :config:`grid.My`

   :Value: 61
   :Option: :opt:`-My`
   :Description: Number of grid points in the y direction.

#. :config:`grid.Mz`

   :Value: 31
   :Option: *no short option*
   :Description: Number of vertical grid levels in the ice.

#. :config:`grid.correct_cell_areas`

   :Value: yes
   :Option: *no short option*
   :Description: Compute corrected cell areas using WGS84 datum (for ice area and volume computations).

#. :config:`grid.ice_vertical_spacing`

   :Value: quadratic
   :Choices: ``quadratic, equal``
   :Option: :opt:`-z_spacing`
   :Description: vertical spacing in the ice

#. :config:`grid.lambda`

   :Value: 4 (pure number)
   :Option: *no short option*
   :Description: Vertical grid spacing parameter. Roughly equal to the factor by which the grid is coarser at an end away from the ice-bedrock interface.

#. :config:`grid.max_stencil_width`

   :Value: 2
   :Option: *no short option*
   :Description: Maximum width of the finite-difference stencil used in PISM.

#. :config:`grid.periodicity`

   :Value: xy
   :Choices: ``none, x, y, xy``
   :Option: :opt:`-periodicity`
   :Description: horizontal grid periodicity

#. :config:`hydrology.cavitation_opening_coefficient`

   :Value: 0.500000 (meter-1)
   :Option: :opt:`-hydrology_cavitation_opening_coefficient`
   :Description: c_1 in notes; coefficient of cavitation opening term in evolution of layer thickness in hydrology::Distributed

#. :config:`hydrology.const_bmelt`

   :Value: 3.168876e-10 (meter / second)
   :Option: :opt:`-hydrology_const_bmelt`
   :Description: default value is equivalent to 1 cm per year of melt; only used if hydrology.use_const_bmelt = 'yes'

#. :config:`hydrology.creep_closure_coefficient`

   :Value: 0.040000 (pure number)
   :Option: :opt:`-hydrology_creep_closure_coefficient`
   :Description: c_2 in notes; coefficient of creep closure term in evolution of layer thickness in hydrology::Distributed

#. :config:`hydrology.gradient_power_in_flux`

   :Value: 1.500000 (pure number)
   :Option: :opt:`-hydrology_gradient_power_in_flux`
   :Description: power `\beta` in Darcy's law `q = - k W^{\alpha} |\nabla \psi|^{\beta-2} \nabla \psi`, for subglacial water layer; used by hydrology::Routing and hydrology::Distributed

#. :config:`hydrology.hydraulic_conductivity`

   :Value: 0.001000 (`m^{2 \beta - \alpha} s^{2 \beta - 3} kg^{1-\beta}`)
   :Option: :opt:`-hydrology_hydraulic_conductivity`
   :Description: = k in notes; lateral conductivity, in Darcy's law, for subglacial water layer; units depend on powers alpha = hydrology.thickness_power_in_flux and beta = hydrology_potential_gradient_power_in_flux; used by hydrology::Routing and hydrology::Distributed

#. :config:`hydrology.maximum_time_step`

   :Value: 1 (years)
   :Option: *no short option*
   :Description: maximum allowed time step length used by hydrology::Routing and hydrology::Distributed

#. :config:`hydrology.model`

   :Value: null
   :Choices: ``null, routing, distributed``
   :Option: :opt:`-hydrology`
   :Description: Basal hydrology sub-model.

#. :config:`hydrology.null_diffuse_till_water`

   :Value: no
   :Option: *no short option*
   :Description: Diffuse stored till water laterally. See equation (11) of :cite:`BBssasliding`

#. :config:`hydrology.null_diffusion_distance`

   :Value: 20000 (meters)
   :Option: *no short option*
   :Description: diffusion distance for till water thickness; see equation (11) in :cite:`BBssasliding`; only active if hydrology.null_diffuse_till_water is set

#. :config:`hydrology.null_diffusion_time`

   :Value: 1000 (years)
   :Option: *no short option*
   :Description: diffusion time for till water thickness; see equation (11) in :cite:`BBssasliding`; only active if hydrology.null_diffuse_till_water is set

#. :config:`hydrology.null_strip_width`

   :Value: -1 (meters)
   :Option: *no short option*
   :Description: if negative then mechanism is inactive; width of strip around computational domain in which water velocity and water amount are set to zero; used by hydrology::Routing and hydrology::Distributed

#. :config:`hydrology.regularizing_porosity`

   :Value: 0.010000 (pure number)
   :Option: :opt:`-hydrology_regularizing_porosity`
   :Description: phi_0 in notes; regularizes pressure equation by multiplying time derivative term

#. :config:`hydrology.roughness_scale`

   :Value: 0.100000 (meters)
   :Option: :opt:`-hydrology_roughness_scale`
   :Description: W_r in notes; roughness scale determining maximum amount of cavitation opening in hydrology::Distributed

#. :config:`hydrology.thickness_power_in_flux`

   :Value: 1.250000 (1)
   :Option: :opt:`-hydrology_thickness_power_in_flux`
   :Description: power `\alpha` in Darcy's law `q = - k W^{\alpha} |\nabla \psi|^{\beta-2} \nabla \psi`, for subglacial water layer; used by hydrology::Routing and hydrology::Distributed

#. :config:`hydrology.tillwat_decay_rate`

   :Value: 3.168876e-11 (meter / second)
   :Option: :opt:`-hydrology_tillwat_decay_rate`
   :Description: default value is equivalent to 1 mm per year; rate at which tillwat is reduced to zero, in absence of other effects like input

#. :config:`hydrology.tillwat_max`

   :Value: 2 (meters)
   :Option: :opt:`-hydrology_tillwat_max`
   :Description: maximum effective thickness of the water stored in till

#. :config:`hydrology.use_const_bmelt`

   :Value: no
   :Option: :opt:`-hydrology_use_const_bmelt`
   :Description: if 'yes', subglacial hydrology model sees basal melt rate which is constant and given by hydrology.const_bmelt

#. :config:`inverse.design.cH1`

   :Value: 0 (1)
   :Option: :opt:`-inv_design_cH1`
   :Description: weight of derivative part of an H1 norm for inversion design variables

#. :config:`inverse.design.cL2`

   :Value: 1 (1)
   :Option: :opt:`-inv_design_cL2`
   :Description: weight of derivative-free part of an H1 norm for inversion design variables

#. :config:`inverse.design.func`

   :Value: sobolevH1
   :Choices: ``sobolevH1, tv``
   :Option: :opt:`-inv_design_func`
   :Description: functional used for inversion design variables

#. :config:`inverse.design.param`

   :Value: exp
   :Choices: ``ident, trunc, square, exp``
   :Option: :opt:`-inv_design_param`
   :Description: parameterization of design variables used during inversion

#. :config:`inverse.design.param_hardav_eps`

   :Value: 10000 (Pascal second^(1/3))
   :Option: *no short option*
   :Description: tiny vertically-averaged hardness used as a substitute for 0 in some tauc parameterizations

#. :config:`inverse.design.param_hardav_scale`

   :Value: 1.000000e+08 (Pascal second^(1/3))
   :Option: *no short option*
   :Description: typical size of ice hardness

#. :config:`inverse.design.param_tauc_eps`

   :Value: 100 (Pascal)
   :Option: *no short option*
   :Description: tiny yield stress used as a substitute for 0 in some tauc parameterizations

#. :config:`inverse.design.param_tauc_scale`

   :Value: 100000 (Pascal)
   :Option: *no short option*
   :Description: typical size of yield stresses

#. :config:`inverse.design.param_trunc_hardav0`

   :Value: 1000000 (Pascal second^(1/3))
   :Option: *no short option*
   :Description: transition point of change to linear behaviour for design variable parameterization type 'trunc'

#. :config:`inverse.design.param_trunc_tauc0`

   :Value: 1000 (Pascal)
   :Option: *no short option*
   :Description: transition point of change to linear behaviour for design variable parameterization type 'trunc'

#. :config:`inverse.log_ratio_scale`

   :Value: 10 (pure number)
   :Option: :opt:`-inv_log_ratio_scale`
   :Description: Reference scale for log-ratio functionals

#. :config:`inverse.ssa.hardav_max`

   :Value: 1.000000e+10 (Pascal second^(1/3))
   :Option: *no short option*
   :Description: Maximum allowed value of hardav for inversions with bound constraints

#. :config:`inverse.ssa.hardav_min`

   :Value: 0 (Pascal second^(1/3))
   :Option: *no short option*
   :Description: Minimum allowed value of hardav for inversions with bound constraints

#. :config:`inverse.ssa.length_scale`

   :Value: 50000 (meters)
   :Option: *no short option*
   :Description: typical length scale for rescaling derivative norms

#. :config:`inverse.ssa.method`

   :Value: tikhonov_lmvm
   :Choices: ``sd, nlcg, ign, tikhonov_lmvm, tikhonov_cg, tikhonov_blmvm, tikhonov_lcl, tikhonov_gn``
   :Option: :opt:`-inv_method`
   :Description: algorithm to use for SSA inversions

#. :config:`inverse.ssa.tauc_max`

   :Value: 5.000000e+07 (Pascal)
   :Option: *no short option*
   :Description: Maximum allowed value of tauc for inversions with bound constraints

#. :config:`inverse.ssa.tauc_min`

   :Value: 0 (Pascal)
   :Option: *no short option*
   :Description: Minimum allowed value of tauc for inversions with bound constraints

#. :config:`inverse.ssa.tv_exponent`

   :Value: 1.200000 (pure number)
   :Option: :opt:`-inv_ssa_tv_exponent`
   :Description: Lebesgue exponent for pseudo-TV norm

#. :config:`inverse.ssa.velocity_eps`

   :Value: 0.100000 (meter / year)
   :Option: *no short option*
   :Description: tiny size of ice velocities during inversion

#. :config:`inverse.ssa.velocity_scale`

   :Value: 100 (meter / year)
   :Option: *no short option*
   :Description: typical size of ice velocities expected during inversion

#. :config:`inverse.state_func`

   :Value: meansquare
   :Choices: ``meansquare, log_ratio, log_relative``
   :Option: :opt:`-inv_state_func`
   :Description: functional used for inversion design variables

#. :config:`inverse.target_misfit`

   :Value: 100 (meter / year)
   :Option: :opt:`-inv_target_misfit`
   :Description: desired root misfit for SSA inversions

#. :config:`inverse.tikhonov.atol`

   :Value: 1.000000e-10 (meter / year)
   :Option: :opt:`-tikhonov_atol`
   :Description: absolute threshold for Tikhonov stopping criterion

#. :config:`inverse.tikhonov.penalty_weight`

   :Value: 1 (1)
   :Option: :opt:`-tikhonov_penalty`
   :Description: penalty parameter for Tikhonov inversion

#. :config:`inverse.tikhonov.ptol`

   :Value: 0.100000 (pure number)
   :Option: :opt:`-tikhonov_ptol`
   :Description: threshold for reaching desired misfit for adaptive Tikhonov algorithms

#. :config:`inverse.tikhonov.rtol`

   :Value: 0.050000 (1)
   :Option: :opt:`-tikhonov_rtol`
   :Description: relative threshold for Tikhonov stopping criterion

#. :config:`ocean.always_grounded`

   :Value: no
   :Option: :opt:`-dry`
   :Description: Dry (ocean-less) simulation; ice is considered grounded regardless of ice thickness, bed elevation, and sea level.

#. :config:`ocean.pik_melt_factor`

   :Value: 0.005000 (1)
   :Option: :opt:`-meltfactor_pik`
   :Description: dimensionless tuning parameter in the '-ocean pik' ocean heat flux parameterization; see :cite:`Martinetal2011`

#. :config:`ocean.runoff_to_ocean_melt_b`

   :Value: 0.150000 (1)
   :Option: *no short option*
   :Description: parameter B in eqn. 1 in :cite:`Aschwanden`

#. :config:`ocean.runoff_to_ocean_melt_power_alpha`

   :Value: 0.540000 (1)
   :Option: *no short option*
   :Description: exponent `\alpha` in eqn. 1 in :cite:`Xu2013`

#. :config:`ocean.runoff_to_ocean_melt_power_beta`

   :Value: 1.170000 (1)
   :Option: *no short option*
   :Description: exponent `\beta` in eqn. 1 in :cite:`Xu2013`

#. :config:`ocean.sub_shelf_heat_flux_into_ice`

   :Value: 0.500000 (W meter-2)
   :Option: *no short option*
   :Description: = J meter-2 second-1; naively chosen default value for heat from ocean; see comments in pism::ocean::Constant::shelf_base_mass_flux().

#. :config:`ocean.three_equation_model_clip_salinity`

   :Value: yes
   :Option: :opt:`-clip_shelf_base_salinity`
   :Description: Clip shelf base salinity so that it is in the range [4, 40] k/kg. See :cite:`HollandJenkins1999`.

#. :config:`output.backup_interval`

   :Value: 1 (hours)
   :Option: :opt:`-backup_interval`
   :Description: wall-clock time between automatic backups

#. :config:`output.backup_size`

   :Value: small
   :Choices: ``none, small, medium, big_2d, big``
   :Option: :opt:`-backup_size`
   :Description: The 'size' of a backup file. See configuration parameters output.sizes.medium, output.sizes.big_2d, output.sizes.big

#. :config:`output.fill_value`

   :Value: -2.000000e+09 (none)
   :Option: *no short option*
   :Description: _FillValue used when saving diagnostic quantities

#. :config:`output.format`

   :Value: netcdf3
   :Choices: ``netcdf3, quilt, netcdf4_parallel, pnetcdf``
   :Option: :opt:`-o_format`
   :Description: The I/O format used for spatial fields; 'netcdf3' is the default, 'netcd4_parallel' is available if PISM was built with parallel NetCDF-4, and 'pnetcdf' is available if PISM was built with PnetCDF.

#. :config:`output.ice_free_thickness_standard`

   :Value: 10 (meters)
   :Option: *no short option*
   :Description: If ice is thinner than this standard then a grid cell is considered ice-free for purposes of reporting glacierized area, volume, etc.

#. :config:`output.runtime.area_scale_factor_log10`

   :Value: 6
   :Option: :opt:`-summary_area_scale_factor_log10`
   :Description: an integer; log base 10 of scale factor to use for area (in km^2) in summary line to stdout

#. :config:`output.file_name`

   :Value: unnamed.nc
   :Option: :opt:`-o`
   :Description: The file to save final model results to.

#. :config:`output.runtime.time_unit_name`

   :Value: year
   :Option: *no short option*
   :Description: Time units used when printing model time, time step, and maximum horizontal velocity at summary to stdout.  Must be valid udunits for time.  (E.g. choose from year,month,day,hour,minute,second.)

#. :config:`output.runtime.time_use_calendar`

   :Value: yes
   :Option: *no short option*
   :Description: Whether to use the current calendar when printing model time in summary to stdout.

#. :config:`output.runtime.viewer.size`

   :Value: 320
   :Option: :opt:`-view_size`
   :Description: default diagnostic viewer size (number of pixels of the longer side)

#. :config:`output.runtime.viewer.variables`

   :Value: *no default*
   :Option: :opt:`-view`
   :Description: comma-separated list of map-plane diagnostic quantities to view at runtime

#. :config:`output.runtime.volume_scale_factor_log10`

   :Value: 6
   :Option: :opt:`-summary_vol_scale_factor_log10`
   :Description: an integer; log base 10 of scale factor to use for volume (in km^3) in summary line to stdout

#. :config:`output.save_size`

   :Value: small
   :Choices: ``none, small, medium, big_2d, big``
   :Option: :opt:`-save_size`
   :Description: The 'size' of a snapshot file. See configuration parameters output.sizes.medium, output.sizes.big_2d, output.sizes.big

#. :config:`output.size`

   :Value: medium
   :Choices: ``none, small, medium, big_2d, big``
   :Option: :opt:`-o_size`
   :Description: The 'size' of an output file. See configuration parameters output.sizes.medium, output.sizes.big_2d, output.sizes.big

#. :config:`output.sizes.big`

   :Value: ``cts, liqfrac, temp, temp_pa, uvel, vvel, wvel, wvel_rel``
   :Option: *no short option*
   :Description: Comma-separated list of variables to write to the output (in addition to model_state variables and variables listed in output.sizes.medium and output.sizes.big_2d) if 'big' output size is selected. Does not include fields written by sub-models.

#. :config:`output.sizes.big_2d`

   :Value: ``age, bfrict, bheatflx, bmelt, bwp, bwprel, cell_area, dbdt, effbwp, enthalpybase, enthalpysurf, flux_divergence, hardav, hydroinput, lat, litho_temp, lon, nuH, ocean_kill_mask, rank, tempbase, tempicethk, tempicethk_basal, temppabase, tempsurf, thk, thksmooth, tillphi, topg, velbar, velbase, wallmelt, wvelbase``
   :Option: *no short option*
   :Description: Comma-separated list of variables to write to the output (in addition to model_state variables and variables listed in output.sizes.medium) if 'big_2d' output size is selected. Does not include fields written by boundary models.

#. :config:`output.sizes.medium`

   :Value: ``bwat, bwatvel, climatic_mass_balance, diffusivity, enthalpy, flux, flux_mag, ice_surface_temp, liqfrac, mask, schoofs_theta, strain_rates, taub_mag, tauc, taud_mag, temp_pa, tillwat, topgsmooth, usurf, velbar_mag, velbase_mag, velsurf, velsurf_mag, wvelsurf``
   :Option: *no short option*
   :Description: Comma-separated list of variables to write to the output (in addition to model_state variables) if 'medium' output size (the default) is selected. Does not include fields written by sub-models.

#. :config:`output.timeseries.buffer_size`

   :Value: 10000
   :Option: *no short option*
   :Description: Number of scalar diagnostic time-series records to hold in memory before writing to disk. (PISM writes this many time-series records to reduce I/O costs.) Send the USR2 signal to flush time-series.

#. :config:`output.timeseries.variables`

   :Value: *no default*
   :Option: :opt:`-ts_vars`
   :Description: Requested scalar (time-series) diagnostics. Leave empty to save all available diagnostics.

#. :config:`output.timeseries.append`

   :Value: false
   :Option: :opt:`-ts_append`
   :Description: If true, append to the scalar time series output file.

#. :config:`output.timeseries.filename`

   :Value: *no default*
   :Option: :opt:`-ts_file`
   :Description: Name of the file to save scalar time series to. Leave empty to disable reporting scalar time-series.

#. :config:`output.variable_order`

   :Value: yxz
   :Choices: ``xyz, yxz, zyx``
   :Option: :opt:`-o_order`
   :Description: Variable order to use in output files.

#. :config:`regional.no_model_strip`

   :Value: 5 (km)
   :Option: :opt:`-no_model_strip`
   :Description: Default width of the 'no model strip' in regional setups.

#. :config:`regional.zero_gradient`

   :Value: false
   :Option: :opt:`-zero_grad_where_no_model`
   :Description: Use zero ice thickness and ice surface gradient in the no_model_mask area.

#. :config:`run_info.institution`

   :Value: *no default*
   :Option: :opt:`-institution`
   :Description: Institution name. This string is written to output files as the 'institution' global attribute.

#. :config:`run_info.title`

   :Value: *no default*
   :Option: :opt:`-title`
   :Description: Free-form string containing a concise description of the current run. This string is written to output files as the 'title' global attribute.

#. :config:`stress_balance.calving_front_stress_bc`

   :Value: no
   :Option: :opt:`-cfbc`
   :Description: Apply CFBC condition as in :cite:`Albrechtetal2011`, :cite:`Winkelmannetal2011`.  May only apply to some stress balances; e.g. SSAFD as of May 2011.  If not set then a strength-extension is used, as in :cite:`BBssasliding`.

#. :config:`stress_balance.ice_free_thickness_standard`

   :Value: 10 (meters)
   :Option: *no short option*
   :Description: If ice is thinner than this standard then a cell is considered ice-free for purposes of computing ice velocity distribution.

#. :config:`stress_balance.model`

   :Value: sia
   :Choices: ``none, prescribed_sliding, sia, ssa, prescribed_sliding+sia, ssa+sia``
   :Option: :opt:`-stress_balance`
   :Description: Stress balance model

#. :config:`stress_balance.sia.Glen_exponent`

   :Value: 3 (pure number)
   :Option: :opt:`-sia_n`
   :Description: Glen exponent in ice flow law for SIA

#. :config:`stress_balance.sia.bed_smoother_range`

   :Value: 5000 (meters)
   :Option: :opt:`-bed_smoother_range`
   :Description: half-width of smoothing domain for stressbalance::BedSmoother, in implementing :cite:`Schoofbasaltopg2003` bed roughness parameterization for SIA; set value to zero to turn off mechanism

#. :config:`stress_balance.sia.e_age_coupling`

   :Value: no
   :Option: :opt:`-e_age_coupling`
   :Description: Couple the SIA enhancement factor to age as in :cite:`Greve`.

#. :config:`stress_balance.sia.enhancement_factor`

   :Value: 1 (1)
   :Option: :opt:`-sia_e`
   :Description: Flow enhancement factor for SIA

#. :config:`stress_balance.sia.enhancement_factor_interglacial`

   :Value: 1 (1)
   :Option: :opt:`-sia_e_interglacial`
   :Description: Flow enhancement factor for SIA; used for ice accumulated during interglacial periods.

#. :config:`stress_balance.sia.flow_law`

   :Value: gpbld
   :Choices: ``arr, arrwarm, gk, gpbld, hooke, isothermal_glen, pb, gpbld3``
   :Option: :opt:`-sia_flow_law`
   :Description: The SIA flow law.

#. :config:`stress_balance.sia.grain_size_age_coupling`

   :Value: no
   :Option: :opt:`-grain_size_age_coupling`
   :Description: Use age of the ice to compute grain size to use with the Goldsby-Kohlstedt :cite:`GoldsbyKohlstedt` flow law

#. :config:`stress_balance.sia.surface_gradient_method`

   :Value: haseloff
   :Choices: ``eta, haseloff, mahaffy``
   :Option: :opt:`-gradient`
   :Description: method used for surface gradient calculation at staggered grid points

#. :config:`stress_balance.ssa.Glen_exponent`

   :Value: 3 (pure number)
   :Option: :opt:`-ssa_n`
   :Description: Glen exponent in ice flow law for SSA

#. :config:`stress_balance.ssa.compute_surface_gradient_inward`

   :Value: no
   :Option: *no short option*
   :Description: If yes then use inward first-order differencing in computing surface gradient in the SSA objects.

#. :config:`stress_balance.ssa.dirichlet_bc`

   :Value: no
   :Option: :opt:`-ssa_dirichlet_bc`
   :Description: apply SSA velocity Dirichlet boundary condition

#. :config:`stress_balance.ssa.enhancement_factor`

   :Value: 1 (1)
   :Option: :opt:`-ssa_e`
   :Description: Flow enhancement factor for SSA

#. :config:`stress_balance.ssa.enhancement_factor_interglacial`

   :Value: 1 (1)
   :Option: :opt:`-ssa_e_interglacial`
   :Description: Flow enhancement factor for SSA; used for ice accumulated during interglacial periods.

#. :config:`stress_balance.ssa.epsilon`

   :Value: 1.000000e+13 (Pascal second meter)
   :Option: :opt:`-ssa_eps`
   :Description: Initial amount of regularization in computation of product of effective viscosity and thickness (`\nu H`).  This default value for `\nu H` comes e.g. from a hardness for the Ross ice shelf (`\bar B`) = 1.9e8 Pa `s^{1/3}` :cite:`MacAyealetal` and a typical strain rate of 0.001 1/year for the Ross ice shelf, giving `\nu = (\bar B) / (2 \cdot 0.001^{2/3})` = 9.49e+14 Pa s ~ 30 MPa year, the value in :cite:`Ritzetal2001`, but with a tiny thickness `H` of about 1 cm.

#. :config:`stress_balance.ssa.fd.brutal_sliding`

   :Value: false
   :Option: :opt:`-brutal_sliding`
   :Description: Enhance sliding speed brutally.

#. :config:`stress_balance.ssa.fd.brutal_sliding_scale`

   :Value: 1 (1)
   :Option: :opt:`-brutal_sliding_scale`
   :Description: Brutal SSA Sliding Scale

#. :config:`stress_balance.ssa.fd.lateral_drag.enabled`

   :Value: false
   :Option: *no short option*
   :Description: set viscosity at ice shelf margin next to ice free bedrock as friction parameterization

#. :config:`stress_balance.ssa.fd.lateral_drag.viscosity`

   :Value: 5.000000e+15 (Pascal second)
   :Option: :opt:`-nu_bedrock`
   :Description: Staggered Viscosity used as side friction parameterization.

#. :config:`stress_balance.ssa.fd.max_iterations`

   :Value: 300
   :Option: :opt:`-ssa_maxi`
   :Description: Maximum number of iterations for the ice viscosity computation, in the SSAFD object

#. :config:`stress_balance.ssa.fd.nuH_iter_failure_underrelaxation`

   :Value: 0.800000 (pure number)
   :Option: :opt:`-ssafd_nuH_iter_failure_underrelaxation`
   :Description: In event of 'Effective viscosity not converged' failure, use outer iteration rule nuH <- nuH + f (nuH - nuH_old), where f is this parameter.

#. :config:`stress_balance.ssa.fd.relative_convergence`

   :Value: 0.000100 (1)
   :Option: :opt:`-ssa_rtol`
   :Description: Relative change tolerance for the effective viscosity in the SSAFD object

#. :config:`stress_balance.ssa.fd.replace_zero_diagonal_entries`

   :Value: yes
   :Option: *no short option*
   :Description: Replace zero diagonal entries in the SSAFD matrix with basal_resistance.beta_ice_free_bedrock to avoid solver failures.

#. :config:`stress_balance.ssa.flow_law`

   :Value: gpbld
   :Choices: ``arr, arrwarm, gpbld, hooke, isothermal_glen, pb, gpbld3``
   :Option: :opt:`-ssa_flow_law`
   :Description: The SSA flow law.

#. :config:`stress_balance.ssa.method`

   :Value: fd
   :Choices: ``fd, fem``
   :Option: :opt:`-ssa_method`
   :Description: Algorithm for computing the SSA solution.

#. :config:`stress_balance.ssa.strength_extension.constant_nu`

   :Value: 9.486807e+14 (Pascal second)
   :Option: *no short option*
   :Description: The SSA is made elliptic by use of a constant value for the product of viscosity (nu) and thickness (H).  This value for nu comes from hardness (bar B)=1.9e8 `Pa s^{1/3}` :cite:`MacAyealetal` and a typical strain rate of 0.001 year-1:  `\nu = (\bar B) / (2 \cdot 0.001^{2/3})`.  Compare the value of 9.45e14 Pa s = 30 MPa year in :cite:`Ritzetal2001`.

#. :config:`stress_balance.ssa.strength_extension.min_thickness`

   :Value: 50 (meters)
   :Option: *no short option*
   :Description: The SSA is made elliptic by use of a constant value for the product of viscosity (nu) and thickness (H).  At ice thicknesses below this value the product nu*H switches from the normal vertical integral to a constant value.  The geometry itself is not affected by this value.

#. :config:`stress_balance.vertical_velocity_approximation`

   :Value: centered
   :Choices: ``centered, upstream``
   :Option: :opt:`-vertical_velocity_approximation`
   :Description: Vertical velocity FD approximation. "Upstream" uses first-order finite difference to compute u_x and v_y. Uses basal velocity to make decisions.

#. :config:`surface.force_to_thickness.alpha`

   :Value: 0.010000 (year-1)
   :Option: *no short option*
   :Description: exponential coefficient in force-to-thickness mechanism

#. :config:`surface.force_to_thickness.ice_free_alpha_factor`

   :Value: 1 (1)
   :Option: *no short option*
   :Description: surface.force_to_thickness.alpha is multiplied by this factor in areas that are ice-free according to the target ice thickness and surface.force_to_thickness.ice_free_thickness_threshold

#. :config:`surface.force_to_thickness.ice_free_thickness_threshold`

   :Value: 1 (meters)
   :Option: *no short option*
   :Description: threshold of ice thickness in the force-to-thickness target field. Used to determine whether to use surface.force_to_thickness.ice_free_alpha_factor.

#. :config:`surface.force_to_thickness.start_time`

   :Value: -4.540000e+09 (years)
   :Option: *no short option*
   :Description: Starting time for the "force to thickness" modifier; the default is "start from the creation of the Earth."

#. :config:`surface.pdd.air_temp_all_precip_as_rain`

   :Value: 275.150000 (Kelvin)
   :Option: *no short option*
   :Description: threshold temperature above which all precipitation is rain; must exceed surface.pdd.air_temp_all_precip_as_snow to avoid division by zero, because difference is in a denominator

#. :config:`surface.pdd.air_temp_all_precip_as_snow`

   :Value: 273.150000 (Kelvin)
   :Option: *no short option*
   :Description: threshold temperature below which all precipitation is snow

#. :config:`surface.pdd.balance_year_start_day`

   :Value: 274
   :Option: *no short option*
   :Description: day of year for October 1st, beginning of the balance year in northern latitudes.

#. :config:`surface.pdd.aschwanden.beta_ice_c`

   :Value: 0.006000 (meter / (Kelvin day))
   :Option: *no short option*
   :Description: water-equivalent thickness; for formula (6) in :cite:`Aschwandenetal2009`

#. :config:`surface.pdd.aschwanden.beta_ice_w`

   :Value: 0.008000 (meter / (Kelvin day))
   :Option: *no short option*
   :Description: water-equivalent thickness; for formula (6) in :cite:`Aschwandenetal2009`

#. :config:`surface.pdd.aschwanden.beta_snow_c`

   :Value: 0.001500 (meter / (Kelvin day))
   :Option: *no short option*
   :Description: water-equivalent thickness; for formula (6) in :cite:`Aschwandenetal2009`

#. :config:`surface.pdd.aschwanden.beta_snow_w`

   :Value: 0.003000 (meter / (Kelvin day))
   :Option: *no short option*
   :Description: water-equivalent thickness; for formula (6) in :cite:`Aschwandenetal2009`

#. :config:`surface.pdd.aschwanden.latitude_beta_w`

   :Value: 77 (degree_north)
   :Option: *no short option*
   :Description: latitude below which to use warm case, in formula (6) in :cite:`Aschwandenetal2009`

#. :config:`surface.pdd.aschwanden.warm_cold_transition_width`

   :Value: 1 (degree_north)
   :Option: *no short option*
   :Description: smoothing width in degrees for linear transition between warm and cold values

#. :config:`surface.pdd.factor_ice`

   :Value: 0.008791 (meter / (Kelvin day))
   :Option: *no short option*
   :Description: EISMINT-Greenland value :cite:`RitzEISMINT`; = (8 mm liquid-water-equivalent) / (pos degree day)

#. :config:`surface.pdd.factor_snow`

   :Value: 0.003297 (meter / (Kelvin day))
   :Option: *no short option*
   :Description: EISMINT-Greenland value :cite:`RitzEISMINT`; = (3 mm liquid-water-equivalent) / (pos degree day)

#. :config:`surface.pdd.fausto.T_c`

   :Value: 272.150000 (Kelvin)
   :Option: *no short option*
   :Description: = -1 + 273.15; for formula (6) in :cite:`Faustoetal2009`

#. :config:`surface.pdd.fausto.T_w`

   :Value: 283.150000 (Kelvin)
   :Option: *no short option*
   :Description: = 10 + 273.15; for formula (6) in :cite:`Faustoetal2009`

#. :config:`surface.pdd.fausto.beta_ice_c`

   :Value: 0.015000 (meter / (Kelvin day))
   :Option: *no short option*
   :Description: water-equivalent thickness; for formula (6) in :cite:`Faustoetal2009`

#. :config:`surface.pdd.fausto.beta_ice_w`

   :Value: 0.007000 (meter / (Kelvin day))
   :Option: *no short option*
   :Description: water-equivalent thickness; for formula (6) in :cite:`Faustoetal2009`

#. :config:`surface.pdd.fausto.beta_snow_c`

   :Value: 0.003000 (meter / (Kelvin day))
   :Option: *no short option*
   :Description: water-equivalent thickness; for formula (6) in :cite:`Faustoetal2009`

#. :config:`surface.pdd.fausto.beta_snow_w`

   :Value: 0.003000 (meter / (Kelvin day))
   :Option: *no short option*
   :Description: water-equivalent thickness; for formula (6) in :cite:`Faustoetal2009`

#. :config:`surface.pdd.fausto.latitude_beta_w`

   :Value: 72 (degree_north)
   :Option: *no short option*
   :Description: latitude below which to use warm case, in formula (6) in :cite:`Faustoetal2009`

#. :config:`surface.pdd.firn_compaction_to_accumulation_ratio`

   :Value: 0.750000 (1)
   :Option: *no short option*
   :Description: How much firn as a fraction of accumulation is turned into ice

#. :config:`surface.pdd.firn_depth_file`

   :Value: *no default*
   :Option: :opt:`-pdd_firn_depth_file`
   :Description: The name of the file to read the firn_depth from.

#. :config:`surface.pdd.temperature_standard_deviation_file`

   :Value: *no default*
   :Option: :opt:`-pdd_sd_file`
   :Description: The name of the file to read air_temp_sd from.

#. :config:`surface.pdd.interpret_precip_as_snow`

   :Value: no
   :Option: *no short option*
   :Description: Interpret precipitation as snow fall.

#. :config:`surface.pdd.max_evals_per_year`

   :Value: 52
   :Option: *no short option*
   :Description: maximum number of times the PDD scheme will ask for air temperature and precipitation to build location-dependent time series for computing (expected) number of positive degree days and snow accumulation; the default means the PDD uses weekly samples of the annual cycle; see also surface.pdd.std_dev

#. :config:`surface.pdd.positive_threshold_temp`

   :Value: 273.150000 (Kelvin)
   :Option: *no short option*
   :Description: temperature used to determine meaning of 'positive' degree day

#. :config:`surface.pdd.refreeze`

   :Value: 0.600000 (1)
   :Option: *no short option*
   :Description: EISMINT-Greenland value :cite:`RitzEISMINT`

#. :config:`surface.pdd.refreeze_ice_melt`

   :Value: yes
   :Option: *no short option*
   :Description: If set to 'yes', refreeze surface.pdd.refreeze fraction of melted ice, otherwise all of the melted ice runs off.

#. :config:`surface.pdd.std_dev`

   :Value: 5 (Kelvin)
   :Option: *no short option*
   :Description: std dev of daily temp variation; = EISMINT-Greenland value :cite:`RitzEISMINT`

#. :config:`surface.pdd.std_dev_lapse_lat_base`

   :Value: 72 (degree_north)
   :Option: *no short option*
   :Description: std_dev is a function of latitude, with value surface.pdd.std_dev at this latitude; this value only active if surface.pdd.std_dev_lapse_lat_rate is nonzero 

#. :config:`surface.pdd.std_dev_lapse_lat_rate`

   :Value: 0 (Kelvin / degree_north)
   :Option: *no short option*
   :Description: std_dev is a function of latitude, with rate of change with respect to latitude given by this constant 

#. :config:`surface.pdd.std_dev_param_a`

   :Value: -0.150000 (pure number)
   :Option: *no short option*
   :Description: Parameter a in Sigma = a*T + b, with T in degrees C. Used only if surface.pdd.std_dev_use_param is set to yes.

#. :config:`surface.pdd.std_dev_param_b`

   :Value: 0.660000 (Kelvin)
   :Option: *no short option*
   :Description: Parameter b in Sigma = a*T + b, with T in degrees C. Used only if surface.pdd.std_dev_use_param is set to yes.

#. :config:`surface.pdd.std_dev_use_param`

   :Value: no
   :Option: *no short option*
   :Description: Parameterize standard deviation as a linear function of air temperature over ice-covered grid cells. The region of application is controlled by geometry.ice_free_thickness_standard.

#. :config:`surface.pressure`

   :Value: 0 (Pascal)
   :Option: *no short option*
   :Description: atmospheric pressure; = pressure at ice surface

#. :config:`surface.temp_to_runoff_a`

   :Value: 0.500000 (K-1)
   :Option: *no short option*
   :Description: a in runoff=a * temp + b

#. :config:`time.calendar`

   :Value: 365_day
   :Choices: ``standard, gregorian, proleptic_gregorian, noleap, 365_day, 360_day, julian, none``
   :Option: :opt:`-calendar`
   :Description: The calendar to use.

#. :config:`time.dimension_name`

   :Value: time
   :Option: *no short option*
   :Description: The name of the time dimension in PISM output files.

#. :config:`time.eemian_end`

   :Value: -114500 (years)
   :Option: *no short option*
   :Description: End of the Eemian interglacial period. See :cite:`Greve97Greenland`.

#. :config:`time.eemian_start`

   :Value: -132000 (years)
   :Option: *no short option*
   :Description: Start of the Eemian interglacial period. See :cite:`Greve97Greenland`.

#. :config:`time.holocene_start`

   :Value: -11000 (years)
   :Option: *no short option*
   :Description: Start of the Holocene interglacial period. See :cite:`Greve97Greenland`.

#. :config:`time.reference_date`

   :Value: 1-1-1
   :Option: *no short option*
   :Description: year-month-day; reference date used for calendar computations and in PISM output files

#. :config:`time.run_length`

   :Value: 1000 (years)
   :Option: *no short option*
   :Description: Default run length

#. :config:`time.start_year`

   :Value: 0 (years)
   :Option: *no short option*
   :Description: Start year.

#. :config:`time_stepping.adaptive_ratio`

   :Value: 0.120000 (1)
   :Option: :opt:`-adapt_ratio`
   :Description: Adaptive time stepping ratio for the explicit scheme for the mass balance equation; :cite:`BBL`, inequality (25)

#. :config:`time_stepping.count_steps`

   :Value: no
   :Option: :opt:`-count_steps`
   :Description: If yes, IceModel::run() will count the number of time steps it took.  Sometimes useful for performance evaluation.  Counts all steps, regardless of whether processes (mass continuity, energy, velocity, ...) occurred within the step.

#. :config:`time_stepping.hit_extra_times`

   :Value: yes
   :Option: :opt:`-extra_force_output_times`
   :Description: Modify the time-stepping mechanism to hit times requested using -extra_times.

#. :config:`time_stepping.hit_multiples`

   :Value: 0 (years)
   :Option: :opt:`-timestep_hit_multiples`
   :Description: Hit every X years, where X is specified using this parameter. Use 0 to disable

#. :config:`time_stepping.hit_save_times`

   :Value: no
   :Option: :opt:`-save_force_output_times`
   :Description: Modify the time-stepping mechanism to hit times requested using -save_times.

#. :config:`time_stepping.hit_ts_times`

   :Value: no
   :Option: *no short option*
   :Description: Modify the time-stepping mechanism to hit times requested using -ts_times.

#. :config:`time_stepping.maximum_time_step`

   :Value: 60 (years)
   :Option: :opt:`-max_dt`
   :Description: Maximum allowed time step length

#. :config:`time_stepping.skip.enabled`

   :Value: no
   :Option: :opt:`-skip`
   :Description: Use the temperature, age, and SSA stress balance computation skipping mechanism.

#. :config:`time_stepping.skip.max`

   :Value: 10
   :Option: :opt:`-skip_max`
   :Description: Number of mass-balance steps, including SIA diffusivity updates, to perform before a the temperature, age, and SSA stress balance computations are done
