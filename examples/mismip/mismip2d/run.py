#!/usr/bin/env python
import MISMIP

# This scripts generates bash scripts that run MISMIP experiments and generates
# all the necessary input files.
#
# Run run.py > my_new_mismip.sh and use that.

try:
    from netCDF4 import Dataset as NC
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

import sys

# The "standard" preamble used in many PISM scripts:
preamble = '''
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
if [ -n "${PISM_BIN:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME     PISM_BIN = $PISM_BIN  (already set)"
else
  PISM_BIN=""    # just a guess
  echo "$SCRIPTNAME     PISM_BIN = $PISM_BIN"
fi
'''


class Experiment:

    "A MISMIP experiment."
    experiment = ""
    mode = 1
    model = 1
    semianalytic = True
    Mx = 151
    My = 3
    Mz = 15
    initials = "ABC"
    executable = "$PISM_DO $PISM_MPIDO $NN ${PISM_BIN}pismr"

    def __init__(self, experiment, model=1, mode=1, Mx=None, Mz=15, semianalytic=True,
                 initials="ABC", executable=None):
        self.model = model
        self.mode = mode
        self.experiment = experiment
        self.initials = initials
        self.semianalytic = semianalytic

        if executable:
            self.executable = executable

        if mode == 3:
            self.Mx = Mx
        else:
            self.Mx = 2 * MISMIP.N(self.mode) + 1

        self.My = 3

        if self.experiment == "2b":
            self.Lz = 7000
        else:
            self.Lz = 6000

    def physics_options(self, input_file, step):
        "Options corresponding to modeling choices."
        config_filename = self.config(step)

        options = ["-energy none",  # isothermal setup; allows selecting cold-mode flow laws
                   "-ssa_flow_law isothermal_glen",  # isothermal setup
                   "-yield_stress constant",
                   "-tauc %e" % MISMIP.C(self.experiment),
                   "-pseudo_plastic",
                   "-gradient eta",
                   "-pseudo_plastic_q %e" % MISMIP.m(self.experiment),
                   "-pseudo_plastic_uthreshold %e" % MISMIP.secpera(),
                   "-calving ocean_kill",  # calving at the present front
                   "-ocean_kill_file %s" % input_file,
                   "-config_override %s" % config_filename,
                   "-ssa_method fd",
                   "-cfbc",                # calving front boundary conditions
                   "-part_grid",           # sub-grid front motion parameterization
                   "-ssafd_ksp_rtol 1e-7",
                   "-ys 0",
                   "-ye %d" % MISMIP.run_length(self.experiment, step),
                   "-options_left",
                   ]

        if self.model == 1:
            options.extend(["-stress_balance ssa"])
        else:
            options.extend(["-stress_balance ssa+sia",
                            "-sia_flow_law isothermal_glen",  # isothermal setup
                            ])

        return options

    def config(self, step):
        '''Generates a config file containing flags and parameters
        for a particular experiment and step.

        This takes care of flags and parameters that *cannot* be set using
        command-line options. (We try to use command-line options whenever we can.)
        '''

        filename = "MISMIP_conf_%s_A%s.nc" % (self.experiment, step)

        nc = NC(filename, 'w', format="NETCDF3_CLASSIC")

        var = nc.createVariable("pism_overrides", 'i')

        attrs = {"ocean.always_grounded": "no",
                 "geometry.update.use_basal_melt_rate": "no",
                 "stress_balance.ssa.compute_surface_gradient_inward": "no",
                 "flow_law.isothermal_Glen.ice_softness": MISMIP.A(self.experiment, step),
                 "constants.ice.density": MISMIP.rho_i(),
                 "constants.sea_water.density": MISMIP.rho_w(),
                 "bootstrapping.defaults.geothermal_flux": 0.0,
                 "stress_balance.ssa.Glen_exponent": MISMIP.n(),
                 "constants.standard_gravity": MISMIP.g(),
                 "ocean.sub_shelf_heat_flux_into_ice": 0.0,
                 }

        if self.model != 1:
            attrs["stress_balance.sia.bed_smoother.range"] = 0.0

        for name, value in attrs.items():
            var.setncattr(name, value)

        nc.close()

        return filename

    def bootstrap_options(self, step):
        boot_filename = "MISMIP_boot_%s_M%s_A%s.nc" % (self.experiment, self.mode, step)

        import prepare
        prepare.pism_bootstrap_file(boot_filename, self.experiment, step, self.mode, N=self.Mx,
                                    semianalytical_profile=self.semianalytic)

        options = ["-i %s -bootstrap" % boot_filename,
                   "-Mx %d" % self.Mx,
                   "-My %d" % self.My,
                   "-Mz %d" % self.Mz,
                   "-Lz %d" % self.Lz]

        return options, boot_filename

    def output_options(self, step):
        output_file = self.output_filename(self.experiment, step)
        extra_file = "ex_" + output_file
        ts_file = "ts_" + output_file

        options = ["-extra_file %s" % extra_file,
                   "-extra_times 0:50:3e4",
                   "-extra_vars thk,topg,velbar_mag,flux_mag,mask,dHdt,usurf,hardav,velbase_mag,nuH,tauc,taud,taub,flux_divergence,cell_grounded_fraction",
                   "-ts_file %s" % ts_file,
                   "-ts_times 0:50:3e4",
                   "-o %s" % output_file,
                   "-o_order zyx",
                   ]

        return output_file, options

    def output_filename(self, experiment, step):
        return "%s%s_%s_M%s_A%s.nc" % (self.initials, self.model, experiment, self.mode, step)

    def options(self, step, input_file=None):
        '''Generates a string of PISM options corresponding to a MISMIP experiment.'''

        if input_file is None:
            input_options, input_file = self.bootstrap_options(step)
        else:
            input_options = ["-i %s" % input_file]

        physics = self.physics_options(input_file, step)

        output_file, output_options = self.output_options(step)

        return output_file, (input_options + physics + output_options)

    def run_step(self, step, input_file=None):
        out, opts = self.options(step, input_file)
        print('echo "# Step %s-%s"' % (self.experiment, step))
        print("%s %s" % (self.executable, ' '.join(opts)))
        print('echo "Done."\n')

        return out

    def run(self, step=None):
        print('echo "# Experiment %s"' % self.experiment)

        if self.experiment in ('1a', '1b'):
            # bootstrap
            input_file = None
            # steps 1 to 9
            steps = list(range(1, 10))

        if self.experiment in ('2a', '2b'):
            # start from step 9 of the corresponding experiment 1
            input_file = self.output_filename(self.experiment.replace("2", "1"), 9)
            # steps 8 to 1
            steps = list(range(8, 0, -1))

        if self.experiment == '3a':
            # bootstrap
            input_file = None
            # steps 1 to 13
            steps = list(range(1, 14))

        if self.experiment == '3b':
            # bootstrap
            input_file = None
            # steps 1 to 15
            steps = list(range(1, 16))

        if step is not None:
            input_file = None
            steps = [step]

        for step in steps:
            input_file = self.run_step(step, input_file)


def run_mismip(initials, executable, semianalytic):
    Mx = 601
    models = (1, 2)
    modes = (1, 2, 3)
    experiments = ('1a', '1b', '2a', '2b', '3a', '3b')

    print(preamble)

    for model in models:
        for mode in modes:
            for experiment in experiments:
                e = Experiment(experiment,
                               initials=initials,
                               executable=executable,
                               model=model, mode=mode, Mx=Mx,
                               semianalytic=semianalytic)
                e.run()


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()

    parser.usage = "%prog [options]"
    parser.description = "Creates a script running MISMIP experiments."
    parser.add_option("--initials", dest="initials", type="string",
                      help="Initials (3 letters)", default="ABC")
    parser.add_option("-e", "--experiment", dest="experiment", type="string",
                      default='1a',
                      help="MISMIP experiments (one of '1a', '1b', '2a', '2b', '3a', '3b')")
    parser.add_option("-s", "--step", dest="step", type="int",
                      help="MISMIP step number")
    parser.add_option("-u", "--uniform_thickness",
                      action="store_false", dest="semianalytic", default=True,
                      help="Use uniform 10 m ice thickness")
    parser.add_option("-a", "--all",
                      action="store_true", dest="all", default=False,
                      help="Run all experiments")
    parser.add_option("-m", "--mode", dest="mode", type="int", default=1,
                      help="MISMIP grid mode")
    parser.add_option("--Mx", dest="Mx", type="int", default=601,
                      help="Custom grid size; use with --mode=3")
    parser.add_option("--model", dest="model", type="int", default=1,
                      help="Models: 1 - SSA only; 2 - SIA+SSA")
    parser.add_option("--executable", dest="executable", type="string",
                      help="Executable to run, e.g. 'mpiexec -n 4 pismr'")

    (opts, args) = parser.parse_args()

    if opts.all:
        run_mismip(opts.initials, opts.executable, opts.semianalytic)
        exit(0)

    def escape(arg):
        if arg.find(" ") >= 0:
            parts = arg.split("=")
            return "%s=\"%s\"" % (parts[0], ' '.join(parts[1:]))
        else:
            return arg

    arg_list = [escape(a) for a in sys.argv]

    print("#!/bin/bash")
    print("# This script was created by examples/mismip/run.py. The command was:")
    print("# %s" % (' '.join(arg_list)))

    if opts.executable is None:
        print(preamble)

    e = Experiment(opts.experiment,
                   initials=opts.initials,
                   executable=opts.executable,
                   model=opts.model,
                   mode=opts.mode,
                   Mx=opts.Mx,
                   semianalytic=opts.semianalytic)

    if opts.step is not None:
        e.run(opts.step)
    else:
        e.run()
