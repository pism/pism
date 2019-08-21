#!/usr/bin/env python

# Copyright (C) 2012, 2013, 2014, 2015 Ricarda Winkelmann, Torsten Albrecht,
# Ed Bueler, and Constantine Khroulev

import numpy as np
import PISMNC
import piktests_utils

# command line arguments
options = piktests_utils.process_options("test_iceberg_removal.nc",
                                         domain_size=1000.0)
p = piktests_utils.Parameters()

# create arrays which will go in output
thk = np.zeros((options.My, options.Mx))  # sheet/shelf thickness
bed = np.zeros_like(thk)                 # bedrock surface elevation
accum = np.zeros_like(thk)                 # accumulation/ ablation
Ts = np.zeros_like(thk) + p.air_temperature

bc_mask = np.zeros_like(thk)
ubar = np.zeros_like(thk)
vbar = np.zeros_like(thk)

dx, dy, x, y = piktests_utils.create_grid(options)


def R(x, y):
    return np.maximum(np.abs(x), np.abs(y))


def C(x, y):
    return np.sqrt(x ** 2 + y ** 2)


if options.square:
    rad = R
else:
    rad = C

# bedrock and ice thickness
for j in range(options.My):
    for i in range(options.Mx):
        radius = rad(x[i], y[j])

        # grounded
        if radius <= p.r_gl and x[i] < 0:
            bed[j, i] = 100.0
            thk[j, i] = p.H0

        # ice shelf
        if radius <= p.r_cf and radius > p.r_gl and options.shelf == True:
            thk[j, i] = (4.0 * p.C / p.Q0 * (radius - p.r_gl) + 1 / p.H0 ** 4) ** (-0.25)

        # ocean (shelf and elsewhere)
        if radius >= p.r_gl or x[i] >= 0:
            accum[j, i] = p.accumulation_rate * p.rho_ice
            bed[j, i] = p.topg_min

# cap ice thickness
thk[thk > p.H0] = p.H0

# set values and locations of Dirichlet boundary conditions
for j in range(options.My):
    for i in range(options.Mx):
        radius = rad(x[i], y[j])
        width = dx * 3

        if radius <= p.r_gl - width and x[i] < 0:
            bc_mask[j, i] = 1.0
        elif radius <= p.r_gl and x[i] < 0:
            bc_mask[j, i] = 1.0
            ubar[j, i] = p.vel_bc * (x[i] / radius)
            vbar[j, i] = p.vel_bc * (y[j] / radius)
        else:
            bc_mask[j, i] = 0.0

ncfile = PISMNC.PISMDataset(options.output_filename, 'w', format='NETCDF3_CLASSIC')
piktests_utils.prepare_output_file(ncfile, x, y)

variables = {"thk": thk,
             "topg": bed,
             "ice_surface_temp": Ts,
             "climatic_mass_balance": accum,
             "bc_mask": bc_mask,
             "u_ssa_bc": ubar,
             "v_ssa_bc": vbar}

piktests_utils.write_data(ncfile, variables)
ncfile.close()

print("Successfully created %s" % options.output_filename)
