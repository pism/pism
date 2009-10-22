#!/usr/bin/env python
from pylab import close, figure, clf, hold, plot, xlabel, ylabel, xticks, yticks, axis, legend, title, grid, show, savefig
from numpy import array, polyfit, polyval, log10, floor, ceil
#from matplotlib.font_manager import FontProperties
import getopt
import sys

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

def plot_errors(nc, x, vars, testname, plot_title, filename = None):

    # This mask lets us choose data corresponding to a particular test:
    test = array(map(chr, nc.variables['test'][:]))
    mask = (test == testname)

    # If we have less than 2 points to plot, then bail.
    if (sum(mask) < 2):
        return

    # Just so that I don't forget that it's log10 and not log.
    f = log10

    # Get the independent variable and transform it. Note that everywhere here
    # I assume that neither dx (dy, dz) not errors can be zero or negative.
    dx = nc.variables[x][mask]
    dim = f(dx)
    
    figure(figsize=(10,6));clf();hold(True)

    colors = ['red', 'blue', 'green', 'black', 'brown', 'cyan']
    for (v,c) in zip(vars,colors):
        # Get a particular variable, transform and fit a line through it:
        data = f(nc.variables[v][mask])
        p = polyfit(dim, data, 1)

        # Try to get the long_name, use short_name if it fails:
        try:
            name = nc.variables[v].long_name
        except:
            name = v

        # Create a label for the independent variable:
        if (x == "dx"):
            dim_name = "\Delta x"
        if (x == "dy"):
            dim_name = "\Delta y"
        if (x == "dz"):
            dim_name = "\Delta z"
        if (x == "dzb"):
            dim_name = "\Delta z_{bed}"

        # Variable label:
        var_label = "%s, $O(%s^{%1.2f})$" % (name, dim_name, p[0])

        # Plot errors and the linear fit:
        plot(dim, data, label=var_label, marker='o', color=c)
        plot(dim, polyval(p, dim), ls="--", color=c)

    # Shrink axes, then expand vertically to have integer powers of 10:
    axis('tight')
    _,_,ymin,ymax = axis()
    axis(ymin = floor(ymin), ymax = ceil(ymax))

    # Round grid spacing in x-ticks:
    xticks(dim, map(lambda(x): "%d" % x, dx))
    # Use default (figured out by matplotlib) locations, but change labels for y-ticks:
    loc,_ = yticks()
    yticks(loc, map(lambda(x) : "$10^{%1.1f}$" % x, loc))

    xlabel("$%s$ (%s)" % (dim_name, nc.variables[x].units))

    # Make sure that all variables given have the same units:
    ylabels = array(map(lambda(x): nc.variables[x].units, vars))
    if (any(ylabels != ylabels[0])):
        print "Incompatible units!"
    else:
        ylabel(ylabels[0])

    # Legend, grid and the title:
    legend(loc='upper left', pad=0.05, labelsep = 0, handletextsep = 0.02, handlelen = 0.02)
    #  prop = FontProperties(size='smaller'),
    grid(True)
    title("Test %s %s (%s)" % (testname, plot_title, nc.source))

    if (filename):
        savefig(filename)

def plot_tests(nc, list_of_tests):
    for test_name in list_of_tests:
        # thickness, volume and eta errors:
        if test_name in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
            plot_errors(nc, 'dx', ["maximum_thickness", "average_thickness"], test_name, "ice thickness errors")
            plot_errors(nc, 'dx', ["volume_percentage"], test_name, "relative ice volume errors")
            plot_errors(nc, 'dx', ["relative_max_eta"], test_name, r"relative max $\eta$ errors")

        # errors that are reported for test E only:
        if (test_name == 'E'):
            plot_errors(nc, 'dx', ["maximum_basal_velocity", "average_basal_velocity"], 'E' , r"basal velocity errors")
            plot_errors(nc, 'dx', ["maximum_ub", "maximum_vb"], 'E' , "basal velocity (ub and vb) errors")
            plot_errors(nc, 'dx', ["relative_basal_velocity"], 'E', "relative basal velocity errors")

        # F and G temperature, sigma and velocity errors:
        if test_name in ['F', 'G']:
            plot_errors(nc, 'dx', ["maximum_sigma", "average_sigma"],
                        test_name, "strain heating errors")
            plot_errors(nc, 'dx', ["maximum_temperature", "average_temperature"],
                        test_name, "ice temperature errors")
            plot_errors(nc, 'dx', ["maximum_basal_temperature", "average_basal_temperature"],
                        test_name, "basal temperature errors")
            plot_errors(nc, 'dx', ["maximum_surface_velocity", "average_surface_velocity"],
                        test_name, "surface velocity errors")
            plot_errors(nc, 'dx', ["maximum_W", "average_W"],
                        test_name, "FIXME velocity errors")

        # I
        # J

        # test K temperature errors:
        if (test_name == 'K'):
            plot_errors(nc, 'dz', ["maximum_temperature", "average_temperature",
                                   "maximum_bedrock_temperature", "average_bedrock_temperature"],
                        'K', "temperature errors")
        # M

def print_help():
    print """Usage:
-i            specifies an input file (which should be a result of running pismv with the -report_file option)
-t,--tests=   specifies a comma-separated list of tests
--help        prints this message
"""

tests_to_plot = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'M']
input = None
try:
    opts, args = getopt.getopt(sys.argv[1:], "i:t:", ["tests=", "help"])
    for opt, arg in opts:
        if opt in ["-t", "--tests"]:
            tests_to_plot = arg.split(',')
        if opt == "-i":
            input = arg
        if opt == "--help":
            print_help()
            sys.exit(0)
    if input:
        nc = NC(input, 'r')
        close('all')
        plot_tests(nc, tests_to_plot)
        show()
    else:
        print_help()
except:
    print "Processing command-line options failed."
    print_help()
    sys.exit(1)

