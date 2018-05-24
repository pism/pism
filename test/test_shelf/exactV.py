#!/usr/bin/env python

from pylab import figure, plot, xlabel, ylabel, title, show, axis, linspace, hold, subplot, grid, step
import numpy as np
import subprocess
import shlex
import sys

from netCDF4 import Dataset as NC

# Setup

secpera = 3.15569259747e7               # seconds per year
rho_sw = 1028.0                         # sea water density
rho_ice = 910.0                         # ice density
standard_gravity = 9.81                 # g
B0 = 1.9e8                              # ice hardness

# "typical constant ice parameter" as defined in the paper and in Van der
# Veen's "Fundamentals of Glacier Dynamics", 1999
C = (rho_ice * standard_gravity * (1.0 - rho_ice / rho_sw) / (4 * B0)) ** 3

# upstream ice thickness
H0 = 600.0                              # meters
# upstream ice velocity
v0 = 300.0 / secpera                      # 300 meters/year
# upstream ice flux
Q0 = H0 * v0


def H(x):
    """Ice thickness."""
    return (4 * C / Q0 * x + 1 / H0 ** 4) ** (-0.25)


def v(x):
    """Ice velocity."""
    return Q0 / H(x)


def x_c(t):
    """Location of the calving front."""
    return Q0 / (4 * C) * ((3 * C * t + 1 / H0 ** 3) ** (4.0 / 3.0) - 1 / H0 ** 4)


def plot_xc(t_years):
    """Plot the location of the calving front."""
    x = x_c(t_years * secpera) / 1000.0   # convert to km
    _, _, y_min, y_max = axis()

    hold(True)
    plot([x, x], [y_min, y_max], '--g')


def run_pismv(Mx, run_length, options, output):
    command = "pismv -test V -y %f -Mx %d %s -o %s" % (run_length, Mx, options, output)
    print("Running %s" % command)
    subprocess.call(shlex.split(command))


def plot_pism_results(filename, figure_title, color, same_figure=False):
    nc = NC(filename)

    time = nc.variables['time'][0] / secpera  # convert to years

    thk = nc.variables['thk'][0, 1, 2:]
    ubar_ssa = nc.variables['velbar_mag'][0, 1, 2:]
    x = nc.variables['x'][:]
    dx = x[1] - x[0]
    Lx = (x[-1] - x[0]) / 2.0
    x_nc = (x[2:] + Lx - 2 * dx) / 1000.0

    hold(True)

    if same_figure == False:
        figure(1)

    subplot(211)
    title(figure_title)
    plotter(x_nc, H(x_nc * 1000.0), color='black', linestyle='dashed')
    plotter(x_nc, thk, color=color, linewidth=2)
    plot_xc(time)
    ylabel("Ice thickness, m")
    axis(xmin=0, xmax=400, ymax=600)
    grid(True)

    subplot(212)
    plotter(x_nc, v(x_nc * 1000.0) * secpera, color='black', linestyle='dashed')
    plotter(x_nc, ubar_ssa, color=color, linewidth=2)
    plot_xc(time)
    axis(xmin=0, xmax=400, ymax=1000)
    xlabel("km")
    ylabel("ice velocity, m/year")
    grid(True)

    nc.close()


import argparse
# Set up the option parser
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.description = """Manages PISM runs reproducing Figure 6 in Albrecht et al
'Parameterization for subgrid-scale motion of ice-shelf calving fronts', 2011"""

parser.epilog = """Model "variants":

0:  no subgrid parameterization, no stress boundary condition at the calving front
1: -cfbc -part_grid
2: -cfbc -part_grid -part_grid_reduce_frontal_thickness

Here -part_grid_reduce_frontal_thickness adjusts the thickness
threshold used to decide when a 'partially filled' cell becomes full.
This is done to try and match the van der Veen profile this particular
setup is based on. Don't use it."""

parser.add_argument("-v", dest="variant", type=int,
                    help="choose the 'model variant', choose from 0, 1, 2", default=2)
parser.add_argument("-Mx", dest="Mx", type=int,
                    help="number of grid points", default=201)
parser.add_argument("-y", dest="y", type=float,
                    help="run length", default=300)
parser.add_argument("-s", dest="step_plot", action='store_true',
                    help="use 'plt.step()' to plot")
options = parser.parse_args()

Mx = options.Mx
x = linspace(0, 400e3, Mx)
run_length = options.y

opt = "-ssa_method fd -Lx 250 -o_order zyx"
extras = " -extra_file ex.nc -extra_vars flux_mag,thk,nuH,flux_divergence,velbar -extra_times 1"

if options.step_plot:
    plotter = step
else:
    plotter = plot

if options.variant == 0:
    run_pismv(Mx, run_length, opt, "out.nc")
    plot_pism_results("out.nc", "Figure 6 (a-b) (control)", 'blue')

    opt = opt + extras
    run_pismv(Mx, run_length, opt + " -max_dt 1", "out.nc")
    plot_pism_results("out.nc", "Figure 6 (a-b) (control)", 'green', same_figure=True)
elif options.variant == 1:
    opt += " -part_grid -cfbc"
    run_pismv(Mx, run_length, opt, "out.nc")
    plot_pism_results("out.nc", "Figure 6 (c-d) (-part_grid)", 'blue')

    opt = opt + extras
    run_pismv(Mx, run_length, opt + " -max_dt 1", "out.nc")
    plot_pism_results("out.nc", "Figure 6 (c-d) (-part_grid)", 'green', same_figure=True)
elif options.variant == 2:
    opt += " -cfbc -part_grid -part_grid_reduce_frontal_thickness"
    run_pismv(Mx, run_length, opt, "out.nc")
    plot_pism_results("out.nc", "Figure 6 (e-f) (-part_grid, reduce frontal thickness)", 'blue')

    opt = opt + extras
    run_pismv(Mx, run_length, opt + " -max_dt 1", "out.nc")
    plot_pism_results("out.nc", "Figure 6 (e-f) (-part_grid)", 'green', same_figure=True)
else:
    print("Wrong variant number. Choose one of 0, 1, 2.")
    sys.exit(1)

show()
