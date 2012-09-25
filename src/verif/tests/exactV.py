#!/usr/bin/env python

# Computes and plots exact solution for "test V", in preparation for
# implementing C version for PISM verification.

from pylab import *
import subprocess

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

def permute(variable, output_order = ('time', 'z', 'zb', 'y', 'x')):
    """Permute dimensions of a NetCDF variable to match the output storage order."""
    input_dimensions = variable.dimensions

    # filter out irrelevant dimensions
    dimensions = filter(lambda(x): x in input_dimensions,
                        output_order)

    # create the mapping
    mapping = map(lambda(x): dimensions.index(x),
                  input_dimensions)

    if mapping:
        return np.transpose(variable[:], mapping)
    else:
        return variable[:]              # so that it does not break processing "mapping"

### Setup

secpera = 3.15569259747e7               # seconds per year
rho_sw = 1028.0                         # sea water density
rho_ice = 910.0                         # ice density
standard_gravity = 9.81                 # g
B0 = 1.9e8                              # ice hardness

# "typical constant ice parameter" as defined in the paper and in Van der
# Veen's "Fundamentals of Glacier Dynamics", 1999
C = (rho_ice * standard_gravity * (1.0 - rho_ice/rho_sw) / (4 * B0))**3

# upstream ice thickness
H0 = 600.0                              # meters
# upstream ice velocity
v0 = 300.0/secpera                      # 300 meters/year
# upstream ice flux
Q0 = H0 * v0;

Mx = 201
x = linspace(0, 400e3, Mx)

def H(x):
    """Ice thickness."""    
    return (4 * C / Q0 * x + 1 / H0**4)**(-0.25)

def v(x):
    """Ice velocity."""    
    return Q0 / H(x)

def x_c(t):
    """Location of the calving front."""
    return Q0 / (4*C) * ((3*C*t + 1/H0**3)**(4.0/3.0) - 1/H0**4)

def plot_xc(t_years):
    """Plot the location of the calving front."""
    x = x_c(t_years * secpera)/1000.0   # convert to km
    a = axis()
    y_min = a[2]
    y_max = a[3]

    hold(True)
    plot([x, x], [y_min, y_max], '--g')

def run_pismv(Mx, run_length, options, output):
    command = "pismv -test V -y %f -Mx %d %s -o %s" % (run_length, Mx, options, output)
    print "Running %s" % command
    subprocess.call(command, shell=True)

def plot_pism_results(figure_number, filename, figure_title, color):
    nc = NC(filename)

    time = nc.variables['time'][0]/secpera # convert to years

    thk = permute(nc.variables['thk'])[0,1,2:]
    ubar_ssa = permute(nc.variables['cbar'])[0,1,2:]
    x = nc.variables['x'][:]
    dx = x[1] - x[0]
    Lx = (x[-1] - x[0]) / 2.0
    x_nc = (x[2:] + Lx - 2*dx) / 1000.0

    hold(True)

    figure(figure_number)

    subplot(211)
    title(figure_title)
    plot(x_nc, H(x_nc*1000.0), color='black', linestyle='dashed')
    plot(x_nc, thk, color=color, linewidth=2)
    plot_xc(time)
    ylabel("Ice thickness, m")
    axis(xmin=0, xmax=400, ymax=600)
    grid(True)

    subplot(212)
    plot(x_nc, v(x_nc*1000.0) * secpera, color='black', linestyle='dashed')
    plot(x_nc, ubar_ssa, color=color, linewidth=2)
    plot_xc(time)
    axis(xmin=0, xmax=400, ymax=1000)
    xlabel("km")
    ylabel("ice velocity, m/year")
    grid(True)

    nc.close()

options = "-ssa_method fd -cfbc -part_grid -Lx 250"

run_pismv(101, 300, options, "out.nc")
plot_pism_results(1, "out.nc", "Figure 6 (b) (-part_grid)", 'blue')

run_pismv(101, 300, options + " -max_dt 1", "out.nc")
plot_pism_results(1, "out.nc", "Figure 6 (b) (-part_grid)", 'green')

run_pismv(101, 300, options + " -part_redist", "out.nc")
plot_pism_results(2, "out.nc", "Figure 6 (c) (-part_grid -part_redist)", 'blue')

run_pismv(101, 300, options + " -part_redist -max_dt 1", "out.nc")
plot_pism_results(2, "out.nc", "Figure 6 (c) (-part_grid -part_redist)", 'green')

show()
