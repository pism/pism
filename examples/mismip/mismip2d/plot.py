#!/usr/bin/env python

import MISMIP

from pylab import figure, subplot, hold, plot, xlabel, ylabel, title, axis, vlines, savefig, text
from sys import exit

import numpy as np
from optparse import OptionParser
import os.path

try:
    from netCDF4 import Dataset as NC
except:
    print("netCDF4 is not installed!")
    sys.exit(1)


def parse_filename(filename, opts):
    "Get MISMIP info from a file name."
    tokens = filename.split('_')
    if tokens[0] == "ex":
        tokens = tokens[1:]

    try:
        model = tokens[0]
        experiment = tokens[1]
        mode = int(tokens[2][1])
        step = int(tokens[3][1])

    except:
        if opts.experiment is None:
            print("Please specify the experiment name (e.g. '-e 1a').")
            exit(0)
        else:
            experiment = opts.experiment

        if opts.step is None:
            print("Please specify the step (e.g. '-s 1').")
            exit(0)
        else:
            step = opts.step

        if opts.model is None:
            print("Please specify the model name (e.g. '-m ABC1').")
            exit(0)
        else:
            model = opts.model

        try:
            nc = NC(filename)
            x = nc.variables['x'][:]
            N = (x.size - 1) / 2
            if N == 150:
                mode = 1
            elif N == 1500:
                mode = 2
            else:
                mode = 3
        except:
            mode = 3

    return model, experiment, mode, step


def process_options():
    "Process command-line options and arguments."
    parser = OptionParser()
    parser.usage = "%prog <input files> [options]"
    parser.description = "Plots the ice flux as a function of the distance from the divide."
    parser.add_option("-o", "--output", dest="output", type="string",
                      help="Output image file name (e.g. -o foo.png)")
    parser.add_option("-e", "--experiment", dest="experiment", type="string",
                      help="MISMIP experiment: 1a,1b,2a,2b,3a,3b (e.g. -e 1a)")
    parser.add_option("-s", "--step", dest="step", type="int",
                      help="MISMIP step: 1,2,3,... (e.g. -s 1)")
    parser.add_option("-m", "--model", dest="model", type="string",
                      help="MISMIP model (e.g. -M ABC1)")
    parser.add_option("-f", "--flux", dest="profile", action="store_false", default=True,
                      help="Plot ice flux only")
    parser.add_option("-p", "--profile", dest="flux", action="store_false", default=True,
                      help="Plot geometry profile only")

    opts, args = parser.parse_args()

    if len(args) == 0:
        print("ERROR: An input file is requied.")
        exit(0)

    if len(args) > 1 and opts.output:
        print("More than one input file given. Ignoring the -o option...\n")
        opts.output = None

    if opts.output and opts.profile and opts.flux:
        print("Please choose between flux (-f) and profile (-p) plots.")
        exit(0)

    return args, opts.output, opts


def read(filename, name):
    "Read a variable and extract the middle row."
    nc = NC(filename)

    try:
        var = nc.variables[name][:]
    except:
        print("ERROR: Variable '%s' not present in '%s'" % (name, filename))
        exit(1)

    N = len(var.shape)
    if N == 1:
        return var[:]               # a coordinate variable ('x')
    elif N == 2:
        return var[1]               # get the middle row
    elif N == 3:
        return var[-1, 1]            # get the middle row of the last record
    else:
        raise Exception("Can't read %s. (It's %d-dimensional.)" % (name, N))


def find_grounding_line(x, topg, thk, mask):
    "Find the modeled grounding line position."
    # "positive" parts of x, topg, thk, mask
    topg = topg[x > 0]
    thk = thk[x > 0]
    mask = mask[x > 0]
    x = x[x > 0]                        # this should go last

    def f(j):
        "See equation (7) in Pattyn et al, 'Role of transition zones in marine ice sheet dynamics', 2005."
        z_sl = 0
        return (z_sl - topg[j]) * MISMIP.rho_w() / (MISMIP.rho_i() * thk[j])

    for j in range(x.size):
        if mask[j] == 2 and mask[j + 1] == 3:  # j is grounded, j+1 floating
            nabla_f = (f(j + 1) - f(j)) / (x[j + 1] - x[j])

            # See equation (8) in Pattyn et al
            return (1.0 - f(j) + nabla_f * x[j]) / nabla_f

    raise Exception("Can't find the grounding line")


def plot_profile(in_file, out_file):
    print("Reading %s to plot geometry profile for model %s, experiment %s, grid mode %s, step %s" % (
        in_file, model, experiment, mode, step))

    if out_file is None:
        out_file = os.path.splitext(in_file)[0] + "-profile.pdf"

    mask = read(in_file, 'mask')
    usurf = read(in_file, 'usurf')
    topg = read(in_file, 'topg')
    thk = read(in_file, 'thk')
    x = read(in_file, 'x')

    # theoretical grounding line position
    xg = MISMIP.x_g(experiment, step)
    # modeled grounding line position
    xg_PISM = find_grounding_line(x, topg, thk, mask)

    # mask out ice-free areas
    usurf = np.ma.array(usurf, mask=mask == 4)

    # compute the lower surface elevation
    lsrf = topg.copy()
    lsrf[mask == 3] = -MISMIP.rho_i() / MISMIP.rho_w() * thk[mask == 3]
    lsrf = np.ma.array(lsrf, mask=mask == 4)

    # convert x to kilometers
    x /= 1e3

    figure(1)
    ax = subplot(111)
    hold(True)
    plot(x, np.zeros_like(x), ls='dotted', color='red')
    plot(x, topg, color='black')
    plot(x, usurf, 'o', color='blue', markersize=4)
    plot(x, lsrf,  'o', color='blue', markersize=4)
    xlabel('distance from the divide, km')
    ylabel('elevation, m')
    title("MISMIP experiment %s, step %d" % (experiment, step))
    text(0.6, 0.9, "$x_g$ (model) = %4.0f km" % (xg_PISM / 1e3), color='r',
         transform=ax.transAxes)
    text(0.6, 0.85, "$x_g$ (theory) = %4.0f km" % (xg / 1e3), color='black',
         transform=ax.transAxes)

    _, _, ymin, ymax = axis(xmin=0, xmax=x.max())
    vlines(xg / 1e3, ymin, ymax, linestyles='dashed', color='black')
    vlines(xg_PISM / 1e3, ymin, ymax, linestyles='dashed', color='red')

    print("Saving '%s'...\n" % out_file)
    savefig(out_file)


def plot_flux(in_file, out_file):
    print("Reading %s to plot ice flux for model %s, experiment %s, grid mode %s, step %s" % (
        in_file, model, experiment, mode, step))

    if out_file is None:
        out_file = os.path.splitext(in_file)[0] + "-flux.pdf"

    x = read(in_file, 'x')
    flux_mag = read(in_file, 'flux_mag')

    # plot positive xs only
    flux_mag = flux_mag[x >= 0]
    x = x[x >= 0]

    figure(2)
    hold(True)

    plot(x / 1e3, flux_mag, 'k.-', markersize=10, linewidth=2)
    plot(x / 1e3, x * MISMIP.a() * MISMIP.secpera(), 'r:', linewidth=1.5)

    title("MISMIP experiment %s, step %d" % (experiment, step))
    xlabel("x ($\mathrm{km}$)", size=14)
    ylabel(r"flux ($\mathrm{m}^2\,\mathrm{a}^{-1}$)", size=14)

    print("Saving '%s'...\n" % out_file)
    savefig(out_file, dpi=300, facecolor='w', edgecolor='w')


if __name__ == "__main__":
    args, out_file, opts = process_options()

    for in_file in args:
        model, experiment, mode, step = parse_filename(in_file, opts)

        if opts.profile:
            plot_profile(in_file, out_file)

        if opts.flux:
            plot_flux(in_file, out_file)
