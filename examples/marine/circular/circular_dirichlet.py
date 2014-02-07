#!/usr/bin/env python

# Copyright (C) 2012, 2013, 2014 Ricarda Winkelmann, Torsten Albrecht,
# Ed Bueler, and Constantine Khroulev

import numpy as np
import PISMNC, piktests_utils

# command line arguments
options = piktests_utils.process_options("circular_dirichlet.nc",
                                         domain_size=1000.0)
p = piktests_utils.Parameters()

dx, dy, x, y = piktests_utils.create_grid(options)

xx, yy = np.meshgrid(x,y)
radius = np.sqrt(xx**2 + yy**2)
# remove the singularity (does not affect the result):
radius[radius == 0] = 1e-16

# Ice thickness
thk = np.zeros((options.My, options.Mx)) # sheet/shelf thickness
if options.shelf:
    thk[radius > p.r_gl] = (4.0 * p.C / p.Q0 * (radius[radius > p.r_gl] - p.r_gl) + 1 / p.H0**4)**(-0.25)
    thk[radius >= p.r_cf] = 0.0
    # cap ice thickness
    thk[thk > p.H0] = p.H0
else:
    thk[radius <= p.r_gl] = p.H0

# Bed topography
bed = np.zeros_like(thk)                 # bedrock surface elevation
bed[radius <= p.r_gl] = -p.H0 * (p.rho_ice / p.rho_ocean) + 1.0
bed[radius >  p.r_gl] = p.topg_min

# Surface mass balance
accum = np.zeros_like(thk)                 # accumulation/ ablation
accum[radius > p.r_gl] = p.accumulation_rate * p.rho_ice # convert to [kg m-2 s-1]

# Surface temperature (irrelevant)
Ts = np.zeros_like(thk) + p.air_temperature

# Dirichlet B.C locations
bcflag = np.zeros_like(thk)
bcflag[radius <= p.r_gl] = 1

# SSA velocity Dirichlet B.C.
ubar = np.zeros_like(thk)
ubar[bcflag == 1] = p.vel_bc * (xx[radius <= p.r_gl] / radius[radius <= p.r_gl])
ubar[bcflag == 0] = 0

vbar = np.zeros_like(thk)
vbar[bcflag == 1] = p.vel_bc * (yy[radius <= p.r_gl] / radius[radius <= p.r_gl])
vbar[bcflag == 0] = 0

ncfile = PISMNC.PISMDataset(options.output_filename, 'w', format='NETCDF3_CLASSIC')
piktests_utils.prepare_output_file(ncfile, x, y)

variables = {"thk"                   : thk,
             "topg"                  : bed,
             "ice_surface_temp"      : Ts,
             "climatic_mass_balance" : accum,
             "bcflag"                : bcflag,
             "u_ssa_bc"              : ubar,
             "v_ssa_bc"              : vbar}

piktests_utils.write_data(ncfile, variables)
ncfile.close()

print "Successfully created %s" % options.output_filename
