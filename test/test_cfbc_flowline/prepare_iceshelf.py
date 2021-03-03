#!/usr/bin/env python

"""This script build a PISM compatible flowline ice shelf setup
   to test the Calving Front Boundary Condition at ice shelves
   and a modified driving stress scheme at the calving front.

   contact: reese@pik-potsdam.de and albrecht@pik-potsdam.de"""


try:
    from netCDF4 import Dataset as NC
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

import numpy as np

# PISM defaults
rhoi=910.0
rhosw=1028
grav=9.81


A0  = 4.6416e-24
#A0  = MISMIP.A('1a', step)
xb  = 100e3

def secpera():
    "Number of seconds per year."
    #return 3.15569259747e7
    return 31536000.0

def L():
    "The length of the MISMIP domain."
    return 1800e3

def bcmask(x):

    bcm = np.zeros_like(x)
    bcm[x < xb] = 1.0

    return np.tile(bcm, (3, 1))


def bed_topography(x):

    bed = np.zeros_like(x)
    #linear
    bed = -(x-xb)*(1000.0)/2000e3 -500.0

    return np.tile(bed, (3, 1))



def surface_mass_balance(x):
    """Computes the surface mass balance."""
    return np.tile(np.zeros_like(x) , (3, 1))


def ice_surface_temp(x):
    """Computes the ice surface temperature (irrelevant)."""
    return np.tile(np.zeros_like(x) + 273.15, (3, 1))


def x(N=None):

    return np.linspace(0, L(), N)


def y(x):
    """Computes y coordinates giving the 1:1 aspect ratio.
    Takes cross-flow grid periodicity into account."""
    dx = x[1] - x[0]
    dy = dx
    return np.array([-dy, 0, dy])


def thickness(x, step, Q0, H0, calving_front=1750e3, pert=False):

    # we expect x to have an odd number of grid points so that one of them is
    # at 0
    if x.size % 2 != 1:
        raise ValueError("x has to have an odd number of points, got %d", x.size)

    # This is basically the analytical solution by van der Veen (1986), 
    # see examples/marine/flowline or test/test_shelf/exactV.py
    dx  = x[1]-x[0]
    thk = np.zeros_like(x)
    B0 = A0**(-1.0/3.0)
    C = (rhoi*grav/(4.0*B0)*(1.0-rhoi/rhosw))**3.0
    thk = (4.0*C*(x-xb)/Q0 + H0**(-4.0))**(-0.25)
    thk[x < xb] = H0
    thk[0]=0.0
    
    if (pert):
      #perutbation i-1
      thk[ x > calving_front-1*dx ] += 50.0 

    thk[x > calving_front] = 0

    return np.tile(thk, (3, 1))



def pism_bootstrap_file(filename, step,
                        v0, H0, calving_front=1750e3, N=None, pert=False):
    import PISMNC

    xx = x(N)
    yy = y(xx)

    print "dx",xx[1]-xx[0]

    v0 = v0/secpera()
    Q0 = H0*v0

    topg = bed_topography(xx)
    thk = thickness(xx, step, Q0, H0, calving_front, pert)
    smb = surface_mass_balance(xx)
    temp = ice_surface_temp(xx)
    vel = Q0/thk*secpera()
    vel[thk==0.0]=0.0
    bcm = bcmask(xx)
    vel0 = np.zeros_like(vel)



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

    nc.define_2d_field('u_ssa_bc',
                       attrs={'units': 'm/yr'
                              })
    nc.write('u_ssa_bc', vel)

    nc.define_2d_field('v_ssa_bc',
                       attrs={'units': 'm/yr'
                              })
    nc.write('v_ssa_bc', vel0)


    nc.define_2d_field('bc_mask',
                       attrs={
                              })
    nc.write('bc_mask', bcm)


    nc.close()


if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    parser.usage = "%prog [options]"
    parser.description = "Creates a flow line iceshelf boostrapping file for use with PISM."
    parser.add_option("-o", dest="output_filename",
                      help="output file")
    parser.add_option("-H", "--thickness", dest="H0", type="float", default=800.0,
                      help="input thickness")
    parser.add_option("-v", "--velocity", dest="v0", type="float", default=300.0,
                      help="input velocity")
    parser.add_option("-s", "--step", dest="step", type="int", default=1,
                      help="MISMIP step")
    parser.add_option("-N", dest="N", type="int", default=3601,
                      help="Custom grid size; use with --mode=3")
    parser.add_option("-c", dest="calving_front", type="float", default=1600e3,
                      help="Calving front location, in meters (e.g. 1600e3)")
    parser.add_option('--perturb_cf', dest='perturb', default=False, action='store_true',
                      help='Add some thickness at last floating front cell')



    (opts, args) = parser.parse_args()


    if not opts.output_filename:
        output_filename = "ICESHELF_%d.nc" % opts.step
    else:
        output_filename = opts.output_filename

    print("Creating ICESHELF setup for step %s in %s..." % (
        opts.step, output_filename))

    pism_bootstrap_file(output_filename,
                        opts.step,
                        opts.v0,
                        opts.H0,
                        calving_front=opts.calving_front,
                        N=opts.N,
                        pert=opts.perturb)

    print("done.")
