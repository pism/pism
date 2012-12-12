#!/usr/bin/env python

## @package vnreport
## Plots error graphs using \c pismv of \c vfnow.py report files.
## 
## Has two command-line options: -i specifies a \c pismv report file to use, \c
## -t selects tests to plot. Use "-t all" to plot error graphs using all the
## error reports saved in a file given with -i.
## 
## For example:
## \verbatim
## vnreport.py -i r1084.nc -t G,B
## \endverbatim
## will use \c r1084.nc and plot tests "G" and "B".

from pylab import close, figure, clf, hold, plot, xlabel, ylabel, xticks, yticks, axis, legend, title, grid, show, savefig
from numpy import array, polyfit, polyval, log10, floor, ceil, unique
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
        print "Skipping Test %s %s (not enough data to plot)" % (testname, plot_title)
        return

    # Get the independent variable and transform it. Note that everywhere here
    # I assume that neither dx (dy, dz) nor errors can be zero or negative.
    dx = nc.variables[x][mask]
    dim = log10(dx)
    
    figure(figsize=(10,6));clf();hold(True)

    colors = ['red', 'blue', 'green', 'black', 'brown', 'cyan']
    for (v,c) in zip(vars,colors):
        # Get a particular variable, transform and fit a line through it:
        data = log10(nc.variables[v][mask])
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

    # Switch to km if dx (dy, dz) are big:
    units = nc.variables[x].units
    if (dx.min() > 1000.0 and (units == "meters")):
        dx = dx / 1000.0
        units = "km"
    # Round grid spacing in x-ticks:
    xticks(dim, map(lambda(x): "%d" % x, dx))
    xlabel("$%s$ (%s)" % (dim_name, units))

    # Use default (figured out by matplotlib) locations, but change labels for y-ticks:
    loc,_ = yticks()
    yticks(loc, map(lambda(x) : "$10^{%1.1f}$" % x, loc))

    # Make sure that all variables given have the same units:
    try:
        ylabels = array(map(lambda(x): nc.variables[x].units, vars))
        if (any(ylabels != ylabels[0])):
            print "Incompatible units!"
        else:
            ylabel(ylabels[0])
    except:
        pass

    # Legend, grid and the title:
    legend(loc='best', borderpad=1, labelspacing = 0.5, handletextpad = 0.75, handlelength = 0.02)
    #  prop = FontProperties(size='smaller'),
    grid(True)
    title("Test %s %s (%s)" % (testname, plot_title, nc.source))

    if (filename):
        savefig(filename)

def plot_tests(nc, list_of_tests):
    for test_name in list_of_tests:
        # thickness, volume and eta errors:
        if test_name in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'L']:
            plot_errors(nc, 'dx', ["maximum_thickness", "average_thickness"], test_name, "ice thickness errors")
            plot_errors(nc, 'dx', ["relative_volume"], test_name, "relative ice volume errors")
            plot_errors(nc, 'dx', ["relative_max_eta"], test_name, r"relative max $\eta$ errors")

        # errors that are reported for test E only:
        if (test_name == 'E'):
            plot_errors(nc, 'dx', ["maximum_basal_velocity", "average_basal_velocity"], 'E' , r"basal velocity errors")
            plot_errors(nc, 'dx', ["maximum_basal_u", "maximum_basal_v"], 'E' , "basal velocity (ub and vb) errors")
            plot_errors(nc, 'dx', ["relative_basal_velocity"], 'E', "relative basal velocity errors")

        # F and G temperature, sigma and velocity errors:
        if test_name in ['F', 'G']:
            plot_errors(nc, 'dx', ["maximum_sigma", "average_sigma"],
                        test_name, "strain heating errors")
            plot_errors(nc, 'dx', ["maximum_temperature", "average_temperature",
                                   "maximum_basal_temperature", "average_basal_temperature"],
                        test_name, "ice temperature errors")

            plot_errors(nc, 'dx', ["maximum_surface_velocity", "maximum_surface_w"],
                        test_name, "maximum ice surface velocity errors")
            plot_errors(nc, 'dx', ["average_surface_velocity", "average_surface_w"],
                        test_name, "average ice surface velocity errors")


        # test I: plot only the u component
        if test_name == 'I':
            plot_errors(nc, 'dy', ["relative_velocity"],
                        test_name, "relative velocity errors")
            plot_errors(nc, 'dy', ["maximum_u", "average_u"],
                        test_name, "velocity errors")
            
        # tests J and M:
        if test_name in ['J', 'M']:
            plot_errors(nc, 'dx', ["relative_velocity"],
                        test_name, "relative velocity errors")
            plot_errors(nc, 'dx', ["maximum_velocity", "maximum_u", "average_u", "maximum_v", "average_v"],
                        test_name, "velocity errors")

        # test K temperature errors:
        if (test_name == 'K'):
            plot_errors(nc, 'dz', ["maximum_temperature", "average_temperature",
                                   "maximum_bedrock_temperature", "average_bedrock_temperature"],
                        'K', "temperature errors")

        # test O temperature and basal melt rate errors:
        if (test_name == 'O'):
            plot_errors(nc, 'dz', ["maximum_temperature", "average_temperature",
                                   "maximum_bedrock_temperature", "average_bedrock_temperature"],
                        'K', "temperature errors")
            plot_errors(nc, 'dz', ["maximum_basal_melt_rate"],
                        'O', "basal melt rate errors")

        # test V: plot only the u component
        if test_name == 'V':
            plot_errors(nc, 'dx', ["relative_velocity"],
                        test_name, "relative velocity errors")
            plot_errors(nc, 'dx', ["maximum_u", "average_u"],
                        test_name, "velocity errors")
def print_help():
    print """Usage:
-i            specifies an input file, a result of running pismv (or vfnow.py)
                with the -report_file option
-t,--tests=   specifies a comma-separated list of tests. Use '-t all'
                to plot all the tests available in the -i file.
--help        prints this message
"""

tests_to_plot = None
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
        available_tests = unique(array(map(chr, nc.variables['test'][:])))

        if len(available_tests) == 1:
            if tests_to_plot == None:
                tests_to_plot = available_tests
        else:
            if (tests_to_plot == None):
                print """Please choose tests to plot using the -t option.
(Input file %s has reports for tests %s available.)""" % (input, str(available_tests))
                sys.exit(0)

        if (tests_to_plot[0] == "all"):
            tests_to_plot = available_tests

        close('all')
        plot_tests(nc, tests_to_plot)
        try:
            # show() will break if we didn't plot anything
            show()
        except:
            pass
    else:
        print_help()
except getopt.GetoptError:
    print "Processing command-line options failed."
    print_help()
    sys.exit(1)

