#!/usr/bin/env python

from pylab import close, figure, clf, hold, plot, xlabel, ylabel, xticks, yticks, axis, legend, title, grid, show, savefig
from numpy import array, polyfit, polyval, log10, floor, ceil, unique
import sys

try:
    from netCDF4 import Dataset as NC
except:
    print("netCDF4 is not installed!")
    sys.exit(1)


class Plotter:

    def __init__(self, save_figures, nc, file_format):
        self.save_figures = save_figures
        self.nc = nc
        self.file_format = file_format

    def plot(self, x, vars, testname, plot_title):
        # This mask lets us choose data corresponding to a particular test:
        test = array(list(map(chr, self.nc.variables['test'][:])))
        mask = (test == testname)

        # If we have less than 2 points to plot, then bail.
        if (sum(mask) < 2):
            print("Skipping Test %s %s (not enough data to plot)" % (testname, plot_title))
            return

        # Get the independent variable and transform it. Note that everywhere here
        # I assume that neither dx (dy, dz) nor errors can be zero or negative.
        dx = self.nc.variables[x][mask]
        dim = log10(dx)

        figure(figsize=(10, 6))
        clf()
        hold(True)

        colors = ['red', 'blue', 'green', 'black', 'brown', 'cyan']
        for (v, c) in zip(vars, colors):
            # Get a particular variable, transform and fit a line through it:
            data = log10(self.nc.variables[v][mask])
            p = polyfit(dim, data, 1)

            # Try to get the long_name, use short_name if it fails:
            try:
                name = self.nc.variables[v].long_name
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

            print("Test {} {}: convergence rate: O(dx^{:1.4f})".format(testname, name, p[0]))

            # Plot errors and the linear fit:
            plot(dim, data, label=var_label, marker='o', color=c)
            plot(dim, polyval(p, dim), ls="--", color=c)

        # Shrink axes, then expand vertically to have integer powers of 10:
        axis('tight')
        _, _, ymin, ymax = axis()
        axis(ymin=floor(ymin), ymax=ceil(ymax))

        # Switch to km if dx (dy, dz) are big:
        units = self.nc.variables[x].units
        if (dx.min() > 1000.0 and (units == "meters")):
            dx = dx / 1000.0
            units = "km"
        # Round grid spacing in x-ticks:
        xticks(dim, ["%d" % x for x in dx])
        xlabel("$%s$ (%s)" % (dim_name, units))

        # Use default (figured out by matplotlib) locations, but change labels for y-ticks:
        loc, _ = yticks()
        yticks(loc, ["$10^{%1.1f}$" % x for x in loc])

        # Make sure that all variables given have the same units:
        try:
            ylabels = array([self.nc.variables[x].units for x in vars])
            if (any(ylabels != ylabels[0])):
                print("Incompatible units!")
            else:
                ylabel(ylabels[0])
        except:
            pass

        # Legend, grid and the title:
        legend(loc='best', borderpad=1, labelspacing=0.5, handletextpad=0.75, handlelength=0.02)
        #  prop = FontProperties(size='smaller'),
        grid(True)
        title("Test %s %s (%s)" % (testname, plot_title, self.nc.source))

        if self.save_figures:
            filename = "%s_%s_%s.%s" % (self.nc.source.replace(" ", "_"),
                                        testname.replace(" ", "_"),
                                        plot_title.replace(" ", "_"),
                                        self.file_format)
            savefig(filename)

    def plot_tests(self, list_of_tests):
        for test_name in list_of_tests:
            # thickness, volume and eta errors:
            if test_name in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'L']:
                self.plot('dx', ["maximum_thickness", "average_thickness"], test_name, "ice thickness errors")
                self.plot('dx', ["relative_volume"], test_name, "relative ice volume errors")
                self.plot('dx', ["relative_max_eta"], test_name, r"relative max eta errors")

            # errors that are reported for test E only:
            if (test_name == 'E'):
                self.plot('dx', ["maximum_basal_velocity", "average_basal_velocity"], 'E', r"basal velocity errors")
                self.plot('dx', ["maximum_basal_u", "maximum_basal_v"], 'E', "basal velocity (ub and vb) errors")
                self.plot('dx', ["relative_basal_velocity"], 'E', "relative basal velocity errors")

            # F and G temperature, sigma and velocity errors:
            if test_name in ['F', 'G']:
                self.plot('dx', ["maximum_sigma", "average_sigma"],
                          test_name, "strain heating errors")
                self.plot('dx', ["maximum_temperature", "average_temperature",
                                 "maximum_basal_temperature", "average_basal_temperature"],
                          test_name, "ice temperature errors")

                self.plot('dx', ["maximum_surface_velocity", "maximum_surface_w"],
                          test_name, "maximum ice surface velocity errors")
                self.plot('dx', ["average_surface_velocity", "average_surface_w"],
                          test_name, "average ice surface velocity errors")

            # test I: plot only the u component
            if test_name == 'I':
                self.plot('dy', ["relative_velocity"],
                          test_name, "relative velocity errors")
                self.plot('dy', ["maximum_u", "average_u"],
                          test_name, "velocity errors")

            # tests J and M:
            if test_name in ['J', 'M']:
                self.plot('dx', ["relative_velocity"],
                          test_name, "relative velocity errors")
                self.plot('dx', ["max_velocity", "maximum_u", "average_u", "maximum_v", "average_v"],
                          test_name, "velocity errors")

            # test K temperature errors:
            if (test_name == 'K'):
                self.plot('dz', ["maximum_temperature", "average_temperature",
                                 "maximum_bedrock_temperature", "average_bedrock_temperature"],
                          'K', "temperature errors")

            # test O temperature and basal melt rate errors:
            if (test_name == 'O'):
                self.plot('dz', ["maximum_temperature", "average_temperature",
                                 "maximum_bedrock_temperature", "average_bedrock_temperature"],
                          'K', "temperature errors")
                self.plot('dz', ["maximum_basal_melt_rate"],
                          'O', "basal melt rate errors")

            # test V: plot only the u component
            if test_name == 'V':
                self.plot('dx', ["relative_velocity"],
                          test_name, "relative velocity errors")
                self.plot('dx', ["maximum_u", "average_u"],
                          test_name, "velocity errors")


from argparse import ArgumentParser
parser = ArgumentParser()
parser.description = """Plot script for PISM verification results."""

parser.add_argument("filename",
                    help="The NetCDF error report file name, usually produces by running vfnow.py")
parser.add_argument("-t", nargs="+", dest="tests_to_plot", default=None,
                    help="Test results to plot (space-delimited list)")
parser.add_argument("--save_figures", dest="save_figures", action="store_true",
                    help="Save figures to .png files")
parser.add_argument("--file_format", dest="file_format", default="png",
                    help="File format for --save_figures (png, pdf, jpg, ...)")

options = parser.parse_args()

input_file = NC(options.filename, 'r')
available_tests = unique(array(list(map(chr, input_file.variables['test'][:]))))
tests_to_plot = options.tests_to_plot

if len(available_tests) == 1:
    if tests_to_plot == None:
        tests_to_plot = available_tests
else:
    if (tests_to_plot == None):
        print("""Please choose tests to plot using the -t option.
(Input file %s has reports for tests %s available.)""" % (input, str(available_tests)))
        sys.exit(0)

if (tests_to_plot[0] == "all"):
    tests_to_plot = available_tests

close('all')

p = Plotter(options.save_figures, input_file, options.file_format)

p.plot_tests(tests_to_plot)
try:
    # show() will break if we didn't plot anything
    if not options.save_figures:
        show()
except:
    pass
