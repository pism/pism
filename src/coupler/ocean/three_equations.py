#!/usr/bin/env python

from numpy import sqrt

shelf_top_surface_temperature = -20.0;  # Celsius
water_latent_heat_fusion      = 3.34e5; # J kg-1
sea_water_density             = 1028.0; # kg m-3
ocean_heat_capacity           = 4170.0; # J kg-1 Kelvin-1
ice_density                   = 910.0;  # kg m-3
ice_heat_capacity             = 2009.0; # J kg-1 Kelvin-1
seconds_per_year              = 3.15569259747e7 # count

def three_equation_model(sea_water_salinity, sea_water_potential_temperature, ice_thickness):
    "Three-equation sub-shelf temperature and melt rate model based on Hellmer et al."

    # Coefficients for linearized freezing point equation for in situ
    # temperature:
    #
    # Tb(salinity, ice_thickness) = a[0] * salinity + a[1] + a[2] * ice_thickness
    a = [-0.0575, 0.0901, -7.61e-4];

    # Coefficients for linearized freezing point equation for potential
    # temperature
    #
    # Theta_b(salinity, ice_thickness) = b[0] * salinity + b[1] + b[2] * ice_thickness
    b = [-0.0575, 0.0921, -7.85e-4];

    # Turbulent heat and salt transfer coefficients:
    gamma_t = 1.00e-4;   # [m/s] RG3417 Default value from Hellmer and Olbers 89
    gamma_s = 5.05e-7;   # [m/s] RG3417 Default value from Hellmer and Olbers 89

    cp_i    = ice_heat_capacity
    cp_w    = ocean_heat_capacity
    L       = water_latent_heat_fusion
    Ts      = shelf_top_surface_temperature
    Sm      = sea_water_salinity
    Theta_m = sea_water_potential_temperature;

    # We solve a quadratic equation for Sb, the salinity at the shelf
    # base.
    #
    # A*Sb^2 + B*Sb + C = 0
    A = a[0] * gamma_s * cp_i - b[0] * gamma_t * cp_w;
    B = (gamma_s * (L - cp_i * (Ts + a[0] * Sm - a[2] * ice_thickness - a[1])) +
         gamma_t * cp_w * (Theta_m - b[2] * ice_thickness - b[1]));
    C = -gamma_s * Sm * (L - cp_i * (Ts - a[2] * ice_thickness - a[1]));
    # Find two roots of the equation:
    S1 = (-B + sqrt(B*B - 4.0 * A * C)) / (2.0 * A);
    S2 = (-B - sqrt(B*B - 4.0 * A * C)) / (2.0 * A);

    if S1 > 0.0:
        basal_salinity = S1
    else:
        basal_salinity = S2;

    shelf_base_temperature = a[0] * basal_salinity + a[1] + a[2] * ice_thickness;

    melt_rate = gamma_s * sea_water_density * (sea_water_salinity - basal_salinity) / (ice_density * basal_salinity);

    return [basal_salinity, shelf_base_temperature, melt_rate*seconds_per_year]
