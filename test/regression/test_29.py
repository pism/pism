#!/usr/bin/env python

import subprocess, shutil, shlex, os
from sys import exit, argv
from netCDF4 import Dataset as NC
import numpy as np

def process_arguments():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("PISM_PATH")
    parser.add_argument("MPIEXEC")
    parser.add_argument("PISM_SOURCE_DIR")

    return parser.parse_args()

def copy_input(opts):
    shutil.copy(os.path.join(opts.PISM_SOURCE_DIR, "test/test_hydrology/inputforP_regression.nc"), ".")

def generate_config():
    """Generates the config file with custom ice softness and hydraulic conductivity."""

    print "generating testPconfig.nc ..."

    nc = NC("testPconfig.nc", 'w')
    pism_overrides = nc.createVariable("pism_overrides", 'b')

    pism_overrides.standard_gravity = 9.81;
    pism_overrides.standard_gravity_doc = "m s-2; = g; acceleration due to gravity on Earth geoid";

    pism_overrides.fresh_water_density = 1000.0;
    pism_overrides.fresh_water_density_doc = "kg m-3; = rhow";

    pism_overrides.ice_softness = 3.1689e-24
    pism_overrides.ice_softness_doc = "Pa-3 s-1; ice softness; NOT DEFAULT"

    pism_overrides.hydrology_hydraulic_conductivity = 1.0e-2 / (1000.0 * 9.81)
    pism_overrides.hydrology_hydraulic_conductivity_doc = "= k; NOT DEFAULT"

    pism_overrides.hydrology_tillwat_max = 0.0;
    pism_overrides.hydrology_tillwat_max_doc = "m; turn off till water mechanism";

    pism_overrides.hydrology_thickness_power_in_flux = 1.0;
    pism_overrides.hydrology_thickness_power_in_flux_doc = "; = alpha in notes";

    pism_overrides.hydrology_potential_gradient_power_in_flux = 2.0;
    pism_overrides.hydrology_potential_gradient_power_in_flux_doc = "; = beta in notes";

    pism_overrides.hydrology_roughness_scale = 1.0;
    pism_overrides.hydrology_roughness_scale_doc = "m; W_r in notes; roughness scale";

    pism_overrides.hydrology_regularizing_porosity = 0.01;
    pism_overrides.hydrology_regularizing_porosity_doc = "[pure]; phi_0 in notes";

    nc.close()

def run_pism(opts):
    cmd = "%s %s/pismr -config_override testPconfig.nc -boot_file inputforP_regression.nc -Mx %d -My %d -Mz 11 -Lz 4000 -hydrology distributed -report_mass_accounting -y 0.08333333333333 -max_dt 0.01 -no_mass -no_energy -ssa_sliding -ssa_dirichlet_bc -o end.nc" % (opts.MPIEXEC, opts.PISM_PATH, 21, 21)

    print cmd
    subprocess.call(shlex.split(cmd))

def check_drift(file1, file2):
    nc1 = NC(file1)
    nc2 = NC(file2)

    stored_drift = {'bwat_max' : 0.024260130035419047,
                    'bwp_max'  : 81469.815066235577,
                    'bwp_avg'  : 6836.2361688114497,
                    'bwat_avg' : 0.0040625953297697154}

    drift = {}
    for name in ("bwat", "bwp"):
        var1 = nc1.variables[name]
        var2 = nc2.variables[name]
        diff = np.abs(np.squeeze(var1[:]) - np.squeeze(var2[:]))

        drift["%s_max" % name] = np.max(diff)
        drift["%s_avg" % name] = np.average(diff)

    print "drift = ", drift
    print "stored_drift = ", stored_drift

    for name in drift.keys():
        rel_diff = np.abs(stored_drift[name] - drift[name]) / stored_drift[name]

        if rel_diff > 1e-3:
            print "Stored and computed drifts in %s differ: %f != %f" % (name, stored_drift[name], drift[name])
            exit(1)

def cleanup():
    for fname in ("inputforP_regression.nc", "testPconfig.nc", "end.nc"):
        os.remove(fname)

if __name__ == "__main__":
    opts = process_arguments()

    print "Copying input files..."
    copy_input(opts)

    print "Generating the -config_override file..."
    generate_config()

    print "Running PISM..."
    run_pism(opts)

    print "Checking the drift..."
    check_drift("inputforP_regression.nc", "end.nc")

    print "Cleaning up..."
    cleanup()
