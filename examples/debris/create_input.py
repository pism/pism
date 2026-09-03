#!/usr/bin/env python3
"""Create the input file for the idealized debris-covered valley glacier of

Y. Verhaegen and P. Huybrechts, "Coupling Debris Transport to 3D Higher-Order Ice Flow
Dynamics to Model the Behavior and Climate Change Response of Debris-Covered Glaciers",
JGR Earth Surface, 2026, doi:10.1029/2025JF008748 (Sections 2.2 and 2.4).

The file contains the bed topography (a headwall merging into a linear bed, with a
U-shaped valley), zero ice thickness, the linear surface mass balance profile, the ice
mask confining the ice to the valley, and the (time-dependent) debris input rate at the
source area near the headwall.
"""

import argparse
import numpy as np
import netCDF4

# bed: h_b(x) = (a - S_b x) + b exp(-x / x_star); valley: - c exp(-(|y - y_c| / w)^m)
a, b, x_star, S_b = 3650.0, 350.0, 150.0, 0.2      # m, m, m, -
c, w, m = 100.0, 250.0, 2.0                        # m, m, -

# mass balance: b_s = beta (x - x_ELA), in m w.e. / yr
beta, x_ELA = -0.003, 1750.0
water_density = 1000.0                             # kg / m^3

# debris source: 2e6 kg / yr over 50 m x 350 m (Table 2), released at t_d = 5 yr
M_in = 2e6                                         # kg / yr
rho_d, phi_d = 2650.0, 0.3
source_length, source_width = 50.0, 350.0          # m
t_d = 5.0                                          # yr

# domain
L_x, L_y = 7500.0, 1500.0                          # m
valley_half_width = 500.0                          # ice mask: |y| < this
headwall_extent = 300.0                            # ice mask: x > this (the source area is ice-free)


def bed_elevation(x, y):
    profile = (a - S_b * x) + b * np.exp(-x / x_star)
    return profile[None, :] - c * np.exp(-(np.abs(y[:, None]) / w)**m)


def create_input(filename, debris_filename, dx):
    x = np.arange(0.5 * dx, L_x, dx)
    y = np.arange(-0.5 * L_y + 0.5 * dx, 0.5 * L_y, dx)
    X, Y = np.meshgrid(x, y)

    topg = bed_elevation(x, y)

    smb = beta * (X - x_ELA) * water_density        # kg m-2 yr-1

    mask = ((np.abs(Y) < valley_half_width) & (X > headwall_extent)).astype(np.float64)

    # debris input: thickness of solid debris per year over the source area, so that
    # the total mass input is M_in
    source = (X > headwall_extent) & (X <= headwall_extent + source_length) & (np.abs(Y) < 0.5 * source_width)
    area = source.sum() * dx * dx
    F_in = M_in / ((1.0 - phi_d) * rho_d * area) if area > 0 else 0.0   # m / yr
    print(f"debris source: {source.sum()} cells, {area:.0f} m^2, F_in = {F_in:.4f} m/yr")

    def coordinates(f):
        f.createDimension("x", len(x))
        f.createDimension("y", len(y))

        v = f.createVariable("x", "f8", ("x",))
        v.units = "m"
        v.axis = "X"
        v.standard_name = "projection_x_coordinate"
        v[:] = x

        v = f.createVariable("y", "f8", ("y",))
        v.units = "m"
        v.axis = "Y"
        v.standard_name = "projection_y_coordinate"
        v[:] = y

    # The bootstrapping file has no time axis: PISM would take the start of the run from
    # the last time record of its input file.
    with netCDF4.Dataset(filename, "w", format="NETCDF4") as f:
        coordinates(f)

        v = f.createVariable("topg", "f8", ("y", "x"))
        v.units = "m"
        v.standard_name = "bedrock_altitude"
        v[:] = topg

        v = f.createVariable("thk", "f8", ("y", "x"))
        v.units = "m"
        v.standard_name = "land_ice_thickness"
        v[:] = 0.0

        v = f.createVariable("climatic_mass_balance", "f8", ("y", "x"))
        v.units = "kg m-2 year-1"
        v.standard_name = "land_ice_surface_specific_mass_balance_flux"
        v[:] = smb

        v = f.createVariable("ice_surface_temp", "f8", ("y", "x"))
        v.units = "kelvin"
        v.long_name = "ice surface temperature (unused: isothermal run)"
        v[:] = 268.15

        v = f.createVariable("land_ice_area_fraction_retreat", "f8", ("y", "x"))
        v.units = "1"
        v.long_name = "maximum ice extent (the valley, excluding the headwall)"
        v[:] = mask

        f.title = "Idealized debris-covered valley glacier (Verhaegen and Huybrechts, 2026)"
        f.source = "examples/debris/create_input.py"

    # the debris input: zero before t_d, F_in on the source area afterwards
    with netCDF4.Dataset(debris_filename, "w", format="NETCDF4") as f:
        coordinates(f)
        f.createDimension("time", None)
        f.createDimension("nv", 2)

        v = f.createVariable("time", "f8", ("time",))
        v.units = "years since 0001-1-1"
        v.calendar = "365_day"
        v.axis = "T"
        v.bounds = "time_bnds"
        v[:] = [0.5 * t_d, t_d + 0.5 * 1e4]

        v = f.createVariable("time_bnds", "f8", ("time", "nv"))
        v[:] = [[0.0, t_d], [t_d, t_d + 1e4]]

        v = f.createVariable("debris_input_rate", "f8", ("time", "y", "x"))
        v.units = "m year-1"
        v.long_name = "rate of debris input from the surrounding terrain (solid debris thickness)"
        v[0, :, :] = 0.0
        v[1, :, :] = np.where(source, F_in, 0.0)

        f.title = "Debris input for the idealized debris-covered valley glacier"
        f.source = "examples/debris/create_input.py"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("output", help="output file name (bootstrapping file)")
    parser.add_argument("--debris-file", default="debris_input.nc",
                        help="name of the debris input (forcing) file (default: debris_input.nc)")
    parser.add_argument("--dx", type=float, default=25.0, help="grid spacing, m (default: 25)")
    args = parser.parse_args()
    create_input(args.output, args.debris_file, args.dx)
