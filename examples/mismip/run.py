#!/usr/bin/env python
import MISMIP

# This scripts generates bash scripts that run MISMIP experiments and generates
# all the necessary input files.
#
# Run run.py > my_new_mismip.sh and use that.

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

# The "standard" preamble used in many PISM scripts:
preamble = """
#!/bin/bash

if [ -n "${SCRIPTNAME:+1}" ] ; then
  echo "[SCRIPTNAME=$SCRIPTNAME (already set)]"
  echo ""
else
  SCRIPTNAME="#(mismip.sh)"
fi

echo
echo "# =================================================================================="
echo "# MISMIP experiments"
echo "# =================================================================================="
echo

set -e  # exit on error

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "mismip.sh 8" then NN = 8
  NN="$1"
fi

echo "$SCRIPTNAME              NN = $NN"

# set MPIDO if using different MPI execution command, for example:
#  $ export PISM_MPIDO="aprun -n "
if [ -n "${PISM_MPIDO:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME      PISM_MPIDO = $PISM_MPIDO  (already set)"
else
  PISM_MPIDO="mpiexec -n "
  echo "$SCRIPTNAME      PISM_MPIDO = $PISM_MPIDO"
fi

# check if env var PISM_DO was set (i.e. PISM_DO=echo for a 'dry' run)
if [ -n "${PISM_DO:+1}" ] ; then  # check if env var DO is already set
  echo "$SCRIPTNAME         PISM_DO = $PISM_DO  (already set)"
else
  PISM_DO=""
fi

# prefix to pism (not to executables)
if [ -n "${PISM_PREFIX:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME     PISM_PREFIX = $PISM_PREFIX  (already set)"
else
  PISM_PREFIX=""    # just a guess
  echo "$SCRIPTNAME     PISM_PREFIX = $PISM_PREFIX"
fi
"""
# " this double-quote makes emacs happier

def physics_options(experiment, step, mode, model=1):
    config_filename = "MISMIP_conf_%s_A%s.nc" % (experiment, step)

    config(config_filename, experiment, step)

    options = ["-cold",                 # allow selecting cold-mode flow laws
               "-flow_law isothermal_glen", # isothermal setup
               "-no_energy",                # isothermal setup
               "-ssa_sliding",              # use SSA
               "-mismip_sliding",           # turn "on" the MISMIP sliding law
               "-ocean_kill",               # calving at the present front
               "-gradient eta", # default method seems to produce artefacts at the grounding line
               "-peridicity y", # periodic in the cross-flow direction
               "-config_override %s" % config_filename,
               "-ssa_method fd",       # use the FD solver that includes the following
               "-cfbc",                # calving front boundary conditions
               "-part_grid",           # sub-grid front motion parameterization
               "-ksp_rtol 1e-7",
               "-config_override %s" % config_filename,
               "-ys 0",
               "-ye 3e4",
               ]

    if model == 1:
        options.extend(["-no_sia"])

    if mode in (2, 3):
        options.extend(["-skip", "-skip_max 10"])

    return options

def config(filename, experiment, step):
    """Generates a config file containing flags and parameters
    for a particular experiment and step.

    This takes care of flags and parameters that *cannot* be set using
    command-line options. (We try to use command-line options whenever we can.)
    """

    nc = NC(filename, 'w')

    var = nc.createVariable("pism_overrides", 'i')

    attrs = {"is_dry_simulation" : "no",
             "include_bmr_in_continuity" : "no",
             "compute_surf_grad_inward_ssa" : "no",
             "surface_gradient_method" : "eta",
             "default_till_phi" : 0.0,
             "ice_softness" : MISMIP.A(experiment, step),
             "ice_density" : MISMIP.rho_i(),
             "sea_water_density" : MISMIP.rho_w(),
             "bootstrapping_geothermal_flux_value_no_var" : 0.0,
             "Glen_exponent" : MISMIP.n(),
             "standard_gravity": MISMIP.g(),
             "ocean_sub_shelf_heat_flux_into_ice" : 0.0,
             "MISMIP_C" : MISMIP.C(experiment),
             "MISMIP_m": MISMIP.m(experiment) }

    for name, value in attrs.iteritems():
        var.setncattr(name, value)

    nc.close()

def bootstrap_options(experiment, step, mode, semianalytical=True, Mz=41):
    Mx = MISMIP.N(mode) + 1
    My = 3

    if experiment == "2b":
        Lz = 7000
    else:
        Lz = 6000

    boot_filename = "MISMIP_boot_%s_M%s_A%s.nc" % (experiment, mode, step)

    import prepare
    prepare.pism_bootstrap_file(boot_filename, experiment, step, mode,
                                semianalytical_profile=semianalytical)

    options = ["-boot_file %s" % boot_filename,
               "-Mx %d" % Mx,
               "-My %d" % My,
               "-Mz %d" % Mz,
               "-Lz %d" % Lz]

    return options


def output_filename(experiment, step, mode, model, initials):
    return "%s%s_%s_M%s_A%s.nc" % (initials, model, experiment, mode, step)

def options(experiment, step, mode, Mz=41, model=1, initials="ABC", semianalytical=True):
    """Generates a string of PISM options corresponding to a MISMIP experiment."""

    input_file = output_filename(experiment, step-1, mode, model, initials)
    output_file = output_filename(experiment, step, mode, model, initials)
    extra_file = "ex_" + output_file

    if step == 1:
        input_options = bootstrap_options(experiment, step, mode, semianalytical, Mz)
    else:
        input_options = ["-i %s" % input_file]

    output_options = ["-extra_file %s" % extra_file,
                      "-extra_times 0:50:3e4",
                      "-extra_vars thk,topg,cbar,velbar,mask",
                      "-o %s" % output_file,
                      "-o_order zyx",
                      ]

    physics = physics_options(experiment, step, mode, model)

    options = input_options + physics + output_options

    return options

def run_experiment(experiment, min_step=1, max_step=15, initials="ABC"):
    print 'echo "# Experiment %s"' % experiment
    for step in range(min_step, max_step + 1):
        print 'echo "# Step %s-%s"' % (experiment, step)

        print ( "$PISM_DO $PISM_MPIDO $NN ${PISM_PREFIX}pismr " +
                ' '.join(options(experiment, step, 2, model=1, initials=initials)) )
        print "echo\n"

def do_mismip(initials):

    print preamble

    for ex in ('1a', '1b'):
        run_experiment(ex, max_step=9, initials=initials)

do_mismip("ABC")
