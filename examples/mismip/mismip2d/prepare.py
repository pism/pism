#!/usr/bin/env python

try:
    from netCDF4 import Dataset as NC
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

import MISMIP

import numpy as np


def bed_topography(experiment, x):
    """Computes bed elevation as a function of x.
    experiment can be '1a', '1b', '2a', '2b', '3a', '3b'.
    """

    return np.tile(-MISMIP.b(experiment, np.abs(x)), (3, 1))


def surface_mass_balance(x):
    """Computes the surface mass balance."""
    return np.tile(np.zeros_like(x) + MISMIP.a(), (3, 1)) * MISMIP.rho_i()


def ice_surface_temp(x):
    """Computes the ice surface temperature (irrelevant)."""
    return np.tile(np.zeros_like(x) + 273.15, (3, 1))


def x(mismip_mode, N=None):
    if mismip_mode in (1, 2):
        return np.linspace(-MISMIP.L(), MISMIP.L(), 2 * MISMIP.N(mismip_mode) + 1)

    return np.linspace(-MISMIP.L(), MISMIP.L(), N)


def y(x):
    """Computes y coordinates giving the 1:1 aspect ratio.
    Takes cross-flow grid periodicity into account."""
    dx = x[1] - x[0]
    dy = dx
    return np.array([-dy, 0, dy])


def thickness(experiment, step, x, calving_front=1750e3, semianalytical_profile=True):

    # we expect x to have an odd number of grid points so that one of them is
    # at 0
    if x.size % 2 != 1:
        raise ValueError("x has to have an odd number of points, got %d", x.size)

    x_nonnegative = x[x >= 0]
    if not semianalytical_profile:
        thk_nonnegative = np.zeros_like(x_nonnegative) + 10
    else:
        thk_nonnegative = MISMIP.thickness(experiment, step, x_nonnegative)

    thk_nonnegative[x_nonnegative > calving_front] = 0

    thk = np.zeros_like(x)
    thk[x >= 0] = thk_nonnegative
    thk[x < 0] = thk_nonnegative[:0:-1]

    return np.tile(thk, (3, 1))


def pism_bootstrap_file(filename, experiment, step, mode,
                        calving_front=1750e3, N=None, semianalytical_profile=True):
    import PISMNC

    xx = x(mode, N)
    yy = y(xx)

    topg = bed_topography(experiment, xx)
    thk = thickness(experiment, step, xx, calving_front, semianalytical_profile)
    smb = surface_mass_balance(xx)
    temp = ice_surface_temp(xx)

    nc = PISMNC.PISMDataset(filename, 'w', format="NETCDF3_CLASSIC")

    nc.create_dimensions(xx, yy)

    nc.define_2d_field('topg',
                       attrs={'units': 'm',
                              'long_name': 'bedrock surface elevation',
                              'standard_name': 'bedrock_altitude'})
    nc.write('topg', topg)

    nc.define_2d_field('thk',
                       attrs={'units': 'm',
                              'long_name': 'ice thickness',
                              'standard_name': 'land_ice_thickness'})
    nc.write('thk', thk)

    nc.define_2d_field('climatic_mass_balance',
                       attrs={'units': 'kg m-2 / s',
                              'long_name': 'ice-equivalent surface mass balance (accumulation/ablation) rate',
                              'standard_name': 'land_ice_surface_specific_mass_balance_flux'})
    nc.write('climatic_mass_balance', smb)

    nc.define_2d_field('ice_surface_temp',
                       attrs={'units': 'Kelvin',
                              'long_name': 'annual average ice surface temperature, below firn processes'})
    nc.write('ice_surface_temp', temp)

    nc.close()


if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    parser.usage = "%prog [options]"
    parser.description = "Creates a MISMIP boostrapping file for use with PISM."
    parser.add_option("-o", dest="output_filename",
                      help="output file")
    parser.add_option("-e", "--experiment", dest="experiment", type="string",
                      help="MISMIP experiment (one of '1a', '1b', '2a', '2b', '3a', '3b')",
                      default="1a")
    parser.add_option("-s", "--step", dest="step", type="int", default=1,
                      help="MISMIP step")
    parser.add_option("-u", "--uniform_thickness",
                      action="store_false", dest="semianalytical_profile", default=True,
                      help="Use uniform 10 m ice thickness")
    parser.add_option("-m", "--mode", dest="mode", type="int", default=2,
                      help="MISMIP grid mode")
    parser.add_option("-N", dest="N", type="int", default=3601,
                      help="Custom grid size; use with --mode=3")
    parser.add_option("-c", dest="calving_front", type="float", default=1600e3,
                      help="Calving front location, in meters (e.g. 1600e3)")

    (opts, args) = parser.parse_args()

    experiments = ('1a', '1b', '2a', '2b', '3a', '3b')
    if opts.experiment not in experiments:
        print("Invalid experiment %s. Has to be one of %s." % (
            opts.experiment, experiments))
        exit(1)

    if not opts.output_filename:
        output_filename = "MISMIP_%s_%d_%d.nc" % (opts.experiment,
                                                  opts.step,
                                                  opts.mode)
    else:
        output_filename = opts.output_filename

    print("Creating MISMIP setup for experiment %s, step %s, grid mode %d in %s..." % (
        opts.experiment, opts.step, opts.mode, output_filename))

    pism_bootstrap_file(output_filename,
                        opts.experiment,
                        opts.step,
                        opts.mode,
                        calving_front=opts.calving_front,
                        N=opts.N,
                        semianalytical_profile=opts.semianalytical_profile)

    print("done.")
