#!/usr/bin/env python

import numpy as np
import argparse

# constants


class Parameters:
    secpera = 3.15569259747e7  # seconds per year
    standard_gravity = 9.81            # m/s^2
    B0 = 1.9e8           # ice hardness

    rho_ice = 910.0               # in kg/m^3
    rho_ocean = 1028.0              # in kg/m^3

    vel_bc = 300     # m/year
    accumulation_rate = 0.3     # m/year
    air_temperature = 247.0   # Kelvin
    domain_size = 1000.0  # km
    topg_min = -3000.0  # m

    # "typical constant ice parameter" as defined in the paper and in Van der
    # Veen's "Fundamentals of Glacier Dynamics", 1999
    C = (rho_ice * standard_gravity * (1.0 - rho_ice / rho_ocean) / (4 * B0)) ** 3
    H0 = 600                     # ice thickness at the grounding line
    Q0 = (vel_bc / secpera) * H0  # flux at the grounding line

    r_gl = 0.25e6             # grounding line radius, meters
    r_cf = 0.45e6             # calving front radius, meters


def process_options(default_filename, domain_size):
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", dest="output_filename",
                        default=default_filename)
    parser.add_argument("-Mx", dest="Mx", type=int, default=301)
    parser.add_argument("-My", dest="My", type=int, default=301)
    parser.add_argument("-domain_size", dest="domain_size", type=float,
                        default=domain_size)
    parser.add_argument("-shelf", dest="shelf", action="store_true")
    parser.add_argument("-square", dest="square", action="store_true")
    return parser.parse_args()


def create_grid(options):
    L = options.domain_size
    dx = L / (options.Mx - 1) * 1000.0  # meters
    dy = L / (options.My - 1) * 1000.0  # meters
    print("dx = %.2f km, dy = %.2f km" % (dx / 1000.0, dy / 1000.0))
    x = np.linspace(-L / 2 * 1000.0, L / 2 * 1000.0, options.Mx)
    y = np.linspace(-L / 2 * 1000.0, L / 2 * 1000.0, options.My)

    return (dx, dy, x, y)


def prepare_output_file(nc, x, y, include_vel_bc=True):
    nc.create_dimensions(x, y)

    attrs = {'long_name': "ice thickness",
             "units": "m",
             "standard_name": "land_ice_thickness",
             "_FillValue": 1.0}
    nc.define_2d_field("thk", attrs=attrs)

    attrs = {'long_name': "bedrock surface elevation",
             "units": "m",
             "standard_name": "bedrock_altitude",
             "_FillValue": -600.0}
    nc.define_2d_field("topg", attrs=attrs)

    attrs = {'long_name': "mean annual temperature at ice surface",
             "units": "Kelvin",
             "_FillValue": 248.0}
    nc.define_2d_field("ice_surface_temp", attrs=attrs)

    ice_density = 910.0
    attrs = {'long_name': "mean annual net ice equivalent accumulation rate",
             "units": "kg m-2 year-1",
             "standard_name": "land_ice_surface_specific_mass_balance_flux",
             "_FillValue": 0.2 * ice_density}  # 0.2 m/year
    nc.define_2d_field("climatic_mass_balance", attrs=attrs)

    if include_vel_bc == False:
        return

    attrs = {'long_name': "Dirichlet boundary condition locations",
             "units": "1",
             "_FillValue": 0}
    nc.define_2d_field("bc_mask", attrs=attrs)

    attrs = {'long_name': "X-component of the SSA velocity boundary conditions",
             "units": "m/year",
             "_FillValue": 0.0}
    nc.define_2d_field("u_ssa_bc", attrs=attrs)

    attrs['long_name'] = "Y-component of the SSA velocity boundary conditions"
    nc.define_2d_field("v_ssa_bc", attrs=attrs)


def write_data(nc, variables):
    for name in list(variables.keys()):
        nc.write(name, variables[name])
