#!/usr/bin/env python
# Copyright (C) 2016-17 Andy Aschwanden

import itertools
from collections import OrderedDict
import numpy as np
import os
import sys
import shlex
from os.path import join, abspath, realpath, dirname

try:
    import subprocess32 as sub
except:
    import subprocess as sub

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys


def current_script_directory():
    import inspect
    filename = inspect.stack(0)[0][1]
    return realpath(dirname(filename))


script_directory = current_script_directory()

from resources import *


def map_dict(val, mdict):
    try:
        return mdict[val]
    except:
        return val


grid_choices = [40, 20, 10]

# set up the option parser
parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.description = "Generating scripts to run Storglaciaren."
parser.add_argument("--restart_file", dest="restart_file",
                    help="Input file to restart from", default=None)
parser.add_argument("-n", '--n_procs', dest="n", type=int,
                    help='''number of cores/processors. default=4.''', default=4)
parser.add_argument("-w", '--wall_time', dest="walltime",
                    help='''walltime. default: 100:00:00.''', default="100:00:00")
parser.add_argument("-q", '--queue', dest="queue", choices=list_queues(),
                    help='''queue. default=debug.''', default='debug')
parser.add_argument("--exstep", dest="exstep",
                    help="Writing interval for spatial time series", default=1)
parser.add_argument("-f", "--o_format", dest="oformat",
                    choices=['netcdf3', 'netcdf4_parallel', 'pnetcdf'],
                    help="output format", default='netcdf3')
parser.add_argument("-g", "--grid", dest="grid", type=int,
                    choices=grid_choices,
                    help="horizontal grid resolution", default=40)
parser.add_argument("--i_dir", dest="input_dir",
                    help="input directory", default=abspath(script_directory))
parser.add_argument("--o_dir", dest="output_dir",
                    help="output directory", default='output')
parser.add_argument("--o_size", dest="osize",
                    choices=['small', 'medium', 'big', 'big_2d', 'custom'],
                    help="output size type", default='custom')
parser.add_argument("-s", "--system", dest="system",
                    choices=list_systems(),
                    help="computer system to use.", default='debug')
parser.add_argument("--spatial_ts", dest="spatial_ts",
                    choices=['default', 'ch', 'none'],
                    help="output size type", default='default')
parser.add_argument("--hydrology", dest="hydrology",
                    choices=['null', 'diffuse', 'routing'],
                    help="Basal hydrology model.", default='diffuse')
parser.add_argument("--stress_balance", dest="stress_balance",
                    choices=['sia', 'ssa+sia', 'ssa'],
                    help="stress balance solver", default='ssa+sia')
parser.add_argument("--vertical_velocity_approximation", dest="vertical_velocity_approximation",
                    choices=['centered', 'upstream'],
                    help="How to approximate vertical velocities", default='upstream')
parser.add_argument("--start_year", dest="start_year", type=int,
                    help="Simulation start year", default=0)
parser.add_argument("--duration", dest="duration", type=int,
                    help="Years to simulate", default=100)
parser.add_argument("--step", dest="step", type=int,
                    help="Step in years for restarting", default=100)
parser.add_argument("-e", "--ensemble_file", dest="ensemble_file",
                    help="File that has all combinations for ensemble study", default='ensemble.csv')
parser.add_argument("-b", "--bootstrap", dest="bootstrap", action="store_true",
                    help="Use bootstrapping", default=False)
parser.add_argument("--climate", dest="climate",
                    choices=['init', 'elev+ftt', 'warming'],
                    help="Climate forcing", default='init')
parser.add_argument("--cryohydrologic_warming", dest="cryohydrologic_warming", action='store_true',
                    help="Turn on cryohydrologic warming", default=False)


options = parser.parse_args()

nn = options.n
input_dir = abspath(options.input_dir)
output_dir = abspath(options.output_dir)

oformat = options.oformat
osize = options.osize
queue = options.queue
walltime = options.walltime
system = options.system

spatial_ts = options.spatial_ts

climate = options.climate
cryohydrologic_warming = options.cryohydrologic_warming
vertical_velocity_approximation = options.vertical_velocity_approximation

exstep = options.exstep
grid = options.grid
hydrology = options.hydrology
stress_balance = options.stress_balance
bootstrap = options.bootstrap

ensemble_file = options.ensemble_file

domain = 'storglaciaren'
pism_exec = generate_domain(domain)

pism_dataname = '$input_dir/pism_storglaciaren_3d.nc'

if options.restart_file is None:
    input_file = pism_dataname
else:
    input_file = options.restart_file

climate_file = pism_dataname

regridvars = 'litho_temp,enthalpy,age,tillwat,bmelt,ice_area_specific_volume,thk'

dirs = {"output": "$output_dir"}
for d in ["performance", "state", "scalar", "spatial", "snap", "jobs", "basins"]:
    dirs[d] = "$output_dir/{dir}".format(dir=d)

if spatial_ts == 'none':
    del dirs["spatial"]

# use the actual path of the run scripts directory (we need it now and
# not during the simulation)
scripts_dir = join(output_dir, "run_scripts")
if not os.path.isdir(scripts_dir):
    os.makedirs(scripts_dir)

# generate the config file *after* creating the output directory
pism_config = 'psg_config'
pism_config_nc = join(output_dir, pism_config + ".nc")

cmd = "ncgen -o {output} {input_dir}/{config}.cdl".format(output=pism_config_nc,
                                                          input_dir=input_dir,
                                                          config=pism_config)
sub.call(shlex.split(cmd))

# these Bash commands are added to the beginning of the run scrips
run_header = """# stop if a variable is not defined
set -u
# stop on errors
set -e

# path to the config file
config="{config}"
# path to the input directory (input data sets are contained in this directory)
input_dir="{input_dir}"
# output directory
output_dir="{output_dir}"

# create required output directories
for each in {dirs};
do
  mkdir -p $each
done

""".format(input_dir=input_dir,
           output_dir=output_dir,
           config=pism_config_nc,
           dirs=" ".join(dirs.values()))

# ########################################################
# set up model initialization
# ########################################################

ssa_n = 3.0
ssa_e = 1.0
tefo = 0.020
phi_min = 5.0
phi_max = 40.
topg_min = -700
topg_max = 700


try:
    combinations = np.loadtxt(ensemble_file, delimiter=',', skiprows=1)
except:
    try:
        combinations = np.genfromtxt(ensemble_file, dtype=None, delimiter=',', skip_header=1, encoding=None)
    except:
        combinations = np.genfromtxt(ensemble_file, dtype=None, delimiter=',', skip_header=1)


tsstep = 'yearly'

scripts = []
scripts_combinded = []

simulation_start_year = options.start_year
simulation_end_year = options.start_year + options.duration
restart_step = options.step

if restart_step > (simulation_end_year - simulation_start_year):
    print('Error:')
    print('restart_step > (simulation_end_year - simulation_start_year): {} > {}'.format(restart_step,
                                                                                         simulation_end_year - simulation_start_year))
    print('Try again')
    import sys
    sys.exit(0)

batch_header, batch_system = make_batch_header(system, nn, walltime, queue)
print(batch_header, batch_system)
for n, combination in enumerate(combinations):

    run_id, sia_e, sia_n, ssa_e, ssa_n, ppq, uth, phi_min, phi_max, topg_min, topg_max = combination

    ttphi = '{},{},{},{}'.format(phi_min, phi_max, topg_min, topg_max)

    name_options = OrderedDict()
    name_options['sb'] = stress_balance
    name_options['climate'] = climate
    try:
        name_options['id'] = '{:03d}'.format(int(run_id))
    except:
        name_options['id'] = '{}'.format(run_id)

    full_exp_name = '_'.join(['_'.join(['_'.join([k, str(v)]) for k, v in name_options.items()])])
    full_outfile = 'g{grid}m_{experiment}.nc'.format(grid=grid, experiment=full_exp_name)

    # All runs in one script file for coarse grids that fit into max walltime
    script_combined = join(scripts_dir, 'sg_g{}m_{}_j.sh'.format(grid, full_exp_name))
    with open(script_combined, 'w') as f_combined:

        outfiles = []
        job_no = 0
        for start in range(simulation_start_year, simulation_end_year, restart_step):
            job_no += 1

            end = start + restart_step

            experiment = '_'.join(['_'.join(['_'.join([k, str(v)])
                                             for k, v in name_options.items()]), '{}'.format(start), '{}'.format(end)])

            script = join(scripts_dir, 'sg_g{}m_{}.sh'.format(grid, experiment))
            scripts.append(script)

            for filename in (script):
                try:
                    os.remove(filename)
                except OSError:
                    pass

            if (start == simulation_start_year):
                f_combined.write(batch_header)
                f_combined.write(run_header)

            with open(script, 'w') as f:

                f.write(batch_header)
                f.write(run_header)

                outfile = '{domain}_g{grid}m_{experiment}.nc'.format(
                    domain=domain.lower(), grid=grid, experiment=experiment)

                pism = generate_prefix_str(pism_exec)

                general_params_dict = {
                    'ys':                          start,
                    'ye':                          end,
                    'calendar':                    '365_day',
                    'climate_forcing_buffer_size': 365,
                    'o':                           join(dirs["state"], outfile),
                    'o_format':                    oformat,
                    'config_override':             "$config"
                }

                if start == simulation_start_year and bootstrap:
                    general_params_dict['bootstrap'] = ''
                    general_params_dict['i'] = pism_dataname
                elif start == simulation_start_year:
                    general_params_dict['bootstrap'] = ''
                    general_params_dict['i'] = pism_dataname
                    general_params_dict['regrid_file'] = input_file
                    general_params_dict['regrid_vars'] = regridvars
                else:
                    general_params_dict['i'] = regridfile

                if osize != 'custom':
                    general_params_dict['o_size'] = osize
                else:
                    general_params_dict['output.sizes.medium'] = 'sftgif,velsurf_mag'

                if start == simulation_start_year:
                    grid_params_dict = generate_grid_description(grid, domain)
                else:
                    grid_params_dict = generate_grid_description(grid, domain, restart=True)

                sb_params_dict = {'sia_e':                              sia_e,
                                  'pseudo_plastic_q':                   ppq,
                                  'pseudo_plastic_threshold':           uth,
                                  'till_effective_fraction_overburden': tefo,
                                  'vertical_velocity_approximation':    vertical_velocity_approximation}

                if start == simulation_start_year:
                    sb_params_dict['topg_to_phi'] = ttphi

                stress_balance_params_dict = generate_stress_balance(stress_balance, sb_params_dict)

                climate_params_dict = generate_climate(
                    climate, **{'energy.ch_warming.enabled': cryohydrologic_warming})

                hydro_params_dict = generate_hydrology(hydrology)

                scalar_ts_dict = generate_scalar_ts(outfile, tsstep,
                                                    start=simulation_start_year,
                                                    end=simulation_end_year,
                                                    odir=dirs["scalar"])

                all_params_dict = merge_dicts(general_params_dict,
                                              grid_params_dict,
                                              stress_balance_params_dict,
                                              climate_params_dict,
                                              hydro_params_dict,
                                              scalar_ts_dict)

                if not spatial_ts == 'none':
                    if spatial_ts == 'default':
                        exvars = default_spatial_ts_vars()
                    else:
                        exvars = ch_spatial_ts_vars()
                    spatial_ts_dict = generate_spatial_ts(outfile, exvars, exstep, odir=dirs["spatial"], split=False)

                    all_params_dict = merge_dicts(all_params_dict,
                                                  spatial_ts_dict)

                all_params = ' \\\n  '.join(["-{} {}".format(k, v) for k, v in all_params_dict.items()])

                if system == 'debug':
                    redirect = ' 2>&1 | tee {jobs}/job_{job_no}'
                else:
                    redirect = ' > {jobs}/job_{job_no}.${job_id} 2>&1'

                template = "{mpido} {pism} {params}" + redirect

                context = merge_dicts(batch_system,
                                      dirs,
                                      {"job_no": job_no, "pism": pism, "params": all_params})
                cmd = template.format(**context)

                f.write(cmd)
                f.write('\n')
                f.write(batch_system.get("footer", ""))

                f_combined.write(cmd)
                f_combined.write('\n\n')

                regridfile = join(dirs["state"], outfile)
                outfiles.append(outfile)

        f_combined.write(batch_system.get("footer", ""))

    scripts_combinded.append(script_combined)


scripts = uniquify_list(scripts)
scripts_combinded = uniquify_list(scripts_combinded)
print('\n'.join([script for script in scripts]))
print('\nwritten\n')
print('\n'.join([script for script in scripts_combinded]))
print('\nwritten\n')
