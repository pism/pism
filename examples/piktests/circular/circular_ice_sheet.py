#!/usr/bin/env python

# Copyright (C) 2012, 2013 Ricarda Winkelmann, Torsten Albrecht,
# Ed Bueler, and Constantine Khroulev

import numpy as np
import scipy.optimize as opt
import PISMNC, piktests_utils

# command line arguments
options = piktests_utils.process_options("circular_withshelf.nc",
                                         domain_size=3600)
p = piktests_utils.Parameters()

dx, dy, x, y = piktests_utils.create_grid(options)

# Compute the calving front radius
x_cf = 1000.0 * options.domain_size / 3 # at 1200 km
y_cf = 1000.0 * options.domain_size / 3 # at 1200 km
r_cf = np.sqrt(x_cf**2 + y_cf**2) # calving front position in m

# create arrays which will go in the output file
thk   = np.zeros((options.My, options.Mx)) # sheet/shelf thickness
bed   = np.zeros_like(thk)                 # bedrock surface elevation
Ts    = np.zeros_like(thk) + p.air_temperature
accum = np.zeros_like(thk) + p.accumulation_rate

def MISMIP_bed(r):
    s  = r / 1000.0
    n  =   729.0
    m2 = -2184.80/(750.0)**2
    m4 =  1031.72/(750.0)**4
    m6 =  -151.72/(750.0)**6
    return n + m2*s**2 + m4*s**4 + m6*s**6 # in m

def MISMIP_thk(r):
    thk_cf  = 200.0             # ice thickness at calving front
    thk_max = 4000.0            # maximal ice thickness in m
    a       = -(thk_cf - thk_max)/(r_cf)**4
    b       = 2*(thk_cf - thk_max)/(r_cf)**2
    c       = thk_max
    return  a * (radius)**4 + b* (radius)**2 + c

# bedrock and ice thickness
for j in range(options.My):
    for i in range(options.Mx):
        radius = np.sqrt(x[i]**2 + y[j]**2) # radius in m

        # set bedrock as in MISMIP experiment
        bed[j,i] = MISMIP_bed(radius)

        # set thickness
        if radius <= r_cf:
            thk[j,i] = MISMIP_thk(radius)

# clip bed topography
bed[bed < p.topg_min] = p.topg_min

# Compute the grounding line radius
I = (options.Mx - 1) / 2
J = (options.My - 1) / 2
H = thk[J,I:]
B = bed[J,I:]
xx = x[I:]
def f(x):
    "floatation criterion: rho_ice/rho_ocean * thk + bed = 0"
    return np.interp(x, xx, p.rho_ice / p.rho_ocean * H + B)

r_gl = opt.bisect(f, xx[0], xx[-1])
print "grounding line radius = %.2f km" % (r_gl/1000.0)

ncfile = PISMNC.PISMDataset(options.output_filename, 'w', format='NETCDF3_CLASSIC')

piktests_utils.prepare_output_file(ncfile, x, y, include_vel_bc=False)

variables = {"thk"                   : thk,
             "topg"                  : bed,
             "ice_surface_temp"      : Ts,
             "climatic_mass_balance" : accum}

piktests_utils.write_data(ncfile, variables)

ncfile.close()
print "Successfully created %s" % options.output_filename
