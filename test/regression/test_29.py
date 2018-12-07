#!/usr/bin/env python

import subprocess
import shutil
import shlex
import os
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

    print("generating testPconfig.nc ...")

    nc = NC("testPconfig.nc", 'w')
    pism_overrides = nc.createVariable("pism_overrides", 'b')

    attrs = {
        "constants.standard_gravity": 9.81,
        "constants.standard_gravity_doc": "m s-2; = g; acceleration due to gravity on Earth geoid",

        "constants.fresh_water.density": 1000.0,
        "constants.fresh_water.density_doc": "kg m-3; = rhow",

        "flow_law.isothermal_Glen.ice_softness": 3.1689e-24,
        "flow_law.isothermal_Glen.ice_softness_doc": "Pa-3 s-1; ice softness; NOT DEFAULT",

        "hydrology.hydraulic_conductivity": 1.0e-2 / (1000.0 * 9.81),
        "hydrology.hydraulic_conductivity_doc": "= k; NOT DEFAULT",

        "hydrology.tillwat_max": 0.0,
        "hydrology.tillwat_max_doc": "m; turn off till water mechanism",

        "hydrology.thickness_power_in_flux": 1.0,
        "hydrology.thickness_power_in_flux_doc": "; = alpha in notes",

        "hydrology.gradient_power_in_flux": 2.0,
        "hydrology.gradient_power_in_flux_doc": "; = beta in notes",

        "hydrology.roughness_scale": 1.0,
        "hydrology.roughness_scale_doc": "m; W_r in notes; roughness scale",

        "hydrology.regularizing_porosity": 0.01,
        "hydrology.regularizing_porosity_doc": "[pure]; phi_0 in notes",

        "basal_yield_stress.model": "constant",
        "basal_yield_stress.model_doc": "only the constant yield stress model works without till",

        "basal_yield_stress.constant.value": 1e6,
        "basal_yield_stress.constant.value_doc": "set default to 'high tauc'"
    }

    for k, v in list(attrs.items()):
        pism_overrides.setncattr(k, v)

    nc.close()


def run_pism(opts):
    cmd = "%s %s/pismr -config_override testPconfig.nc -i inputforP_regression.nc -bootstrap -Mx %d -My %d -Mz 11 -Lz 4000 -hydrology distributed -report_mass_accounting -y 0.08333333333333 -max_dt 0.01 -no_mass -energy none -stress_balance ssa+sia -ssa_dirichlet_bc -o end.nc" % (
        opts.MPIEXEC, opts.PISM_PATH, 21, 21)

    print(cmd)
    subprocess.call(shlex.split(cmd))


def check_drift(file1, file2):
    nc1 = NC(file1)
    nc2 = NC(file2)

    stored_drift = {'bwat_max': 0.023524626576411189,
                    'bwp_max':  79552.478734239354,
                    'bwp_avg':  6261.1642337484445,
                    'bwat_avg': 0.0034449380393343091}

    drift = {}
    for name in ("bwat", "bwp"):
        var1 = nc1.variables[name]
        var2 = nc2.variables[name]
        diff = np.abs(np.squeeze(var1[:]) - np.squeeze(var2[:]))

        drift["%s_max" % name] = np.max(diff)
        drift["%s_avg" % name] = np.average(diff)

    print("drift        = ", drift)
    print("stored_drift = ", stored_drift)

    for name in list(drift.keys()):
        rel_diff = np.abs(stored_drift[name] - drift[name]) / stored_drift[name]

        if rel_diff > 1e-3:
            print("Stored and computed drifts in %s differ: %f != %f" % (name, stored_drift[name], drift[name]))
            exit(1)


def cleanup():
    for fname in ("inputforP_regression.nc", "testPconfig.nc", "end.nc"):
        os.remove(fname)


if __name__ == "__main__":
    opts = process_arguments()

    print("Copying input files...")
    copy_input(opts)

    print("Generating the -config_override file...")
    generate_config()

    print("Running PISM...")
    run_pism(opts)

    print("Checking the drift...")
    check_drift("inputforP_regression.nc", "end.nc")

    print("Cleaning up...")
    cleanup()
