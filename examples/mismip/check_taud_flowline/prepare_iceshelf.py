#!/usr/bin/env python3

from netCDF4 import Dataset as NC

import MISMIP

import numpy as np

# the constant used to protect from division by zero and raising negative numbers to a
# negative power
eps = 1e-16

def bcmask(x):

    bcm = np.zeros_like(x)

    bcm[x < 100e3] = 1.0
    #bcm[x < 105e3] = 1.0

    return np.tile(bcm, (3, 1))


def bed_topography(x):

    bed = np.zeros_like(x)

    #linear
    bed = -(x-100e3)*(1000.0)/2000e3 -500.0
    #constant
    #bed[x >= 100e3] = -500.0
    #bed[x < 100e3] = -500.0


    return np.tile(bed, (3, 1))


def surface_mass_balance(x):
    """Computes the surface mass balance."""
    return np.tile(np.zeros_like(x) + MISMIP.a()*0.0, (3, 1)) * MISMIP.rho_i()


def ice_surface_temp(x):
    """Computes the ice surface temperature (irrelevant)."""
    return np.tile(np.zeros_like(x) + 273.15, (3, 1))


def x(N=None):
    return np.linspace(0, MISMIP.L(), N)


def y(x):
    """Computes y coordinates giving the 1:1 aspect ratio.
    Takes cross-flow grid periodicity into account."""
    dx = x[1] - x[0]
    dy = dx
    return np.array([-dy, 0, dy])


def thickness(x, step, Q0, H0, calving_front=1750e3, perturbation='default'):

    # we expect x to have an odd number of grid points so that one of them is
    # at 0
    if x.size % 2 != 1:
        raise ValueError("x has to have an odd number of points, got %d", x.size)

    dx=x[1]-x[0]

    A0 = MISMIP.A('1a', step)
    print(A0,A0**(-1.0/3.0))
    B0 = A0**(-1.0/MISMIP.n())
    #C = (900.0*9.8/(4.0*B0)*(1.0-900.0/1000.0))**3.0
    C = (MISMIP.rho_i()*MISMIP.g()/(4.0*B0)*(1.0-MISMIP.rho_i()/MISMIP.rho_w()))**MISMIP.n()

    thk = (4.0*C*np.maximum(x-100e3, eps)/Q0 + H0**(-4.0))**(-0.25)
    #thk = (x-100e3)*(-1e-4) + 800.0
    thk[x < 100e3] = 800.0
    thk[0]=0.0
    if perturbation=='p1':
      #perturbation i
      thk[ x > calving_front-1*dx ] -= 50.0
    elif perturbation=='p2':
      #perutbation i-1
      thk[ x > calving_front-2*dx ] -= 50.0 
      thk[ x > calving_front-1*dx ] += 50.0
    elif perturbation=='p3':
      #perutbation i-2
      thk[ x > calving_front-3*dx ] -= 50.0
      thk[ x > calving_front-2*dx ] += 50.0 
    #perutbation gli+1
    #thk[ x > 100e3 + 0.5*dx ] -= 200.0
    #thk[ x > 100e3+1.5*dx ] += 200.0
    #perutbation gli+2
    #thk[ x > 100e3 + 1.5*dx ] -= 200.0
    #thk[ x > 100e3+2.5*dx ] += 200.0

    thk[x > calving_front] = 0.0

    return np.tile(thk, (3, 1))

def ocean_kill(H):
    """Defines a ice front boundary, at which ice calves off."""
    ok=np.ones_like(H)
    ok[H==0]=0
    return ok


def config(step):
        '''Generates a config file containing flags and parameters
        for a particular experiment and step.

        This takes care of flags and parameters that *cannot* be set using
        command-line options. (We try to use command-line options whenever we can.)
        '''
        model = 1
        experiment = '1a'
        filename = "MISMIP_conf_%s_A%s.nc" % (experiment, step)

        nc = NC(filename, 'w', format="NETCDF3_CLASSIC")

        var = nc.createVariable("pism_overrides", 'i')

        attrs = {"stress_balance.ssa.fd.flow_line_mode": "true",
                 "geometry.update.use_basal_melt_rate": "no",
                 "stress_balance.ssa.compute_surface_gradient_inward": "no",
                 "flow_law.isothermal_Glen.ice_softness": MISMIP.A(experiment, step),
                 "constants.ice.density": MISMIP.rho_i(),
                 "constants.sea_water.density": MISMIP.rho_w(),
                 "bootstrapping.defaults.geothermal_flux": 0.0,
                 "stress_balance.ssa.Glen_exponent": MISMIP.n(),
                 "constants.standard_gravity": MISMIP.g(),
                 "ocean.sub_shelf_heat_flux_into_ice": 0.0,
                 }

        if model != 1:
            attrs["stress_balance.sia.bed_smoother.range"] = 0.0

        for name, value in attrs.items():
            var.setncattr(name, value)

        nc.close()

        return filename



def pism_bootstrap_file(filename, step,
                        v0, H0, calving_front=1750e3, N=None,p="default"):
    import PISMNC

    xx = x(N)
    yy = y(xx)

    print("dx",xx[1]-xx[0])

    v0 = v0/MISMIP.secpera()
    #H0 = 800.0
    Q0 = H0*v0

    topg = bed_topography(xx)
    thk = thickness(xx, step, Q0, H0, calving_front, p)
    smb = surface_mass_balance(xx)
    temp = ice_surface_temp(xx)

    vel = Q0/np.maximum(thk, eps)*MISMIP.secpera()
    vel[thk==0.0]=0.0

    #vel = np.nan_to_num(Q0/thk)*MISMIP.secpera()
    bcm = bcmask(xx)
    vel0 = np.zeros_like(vel)
    okill = ocean_kill(thk)





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
                       attrs={'units': 'kg m^-2 / s',
                              'long_name': 'ice-equivalent surface mass balance (accumulation/ablation) rate',
                              'standard_name': 'land_ice_surface_specific_mass_balance_flux'})
    nc.write('climatic_mass_balance', smb)

    nc.define_2d_field('ice_surface_temp',
                       attrs={'units': 'kelvin',
                              'long_name': 'annual average ice surface temperature, below firn processes'})
    nc.write('ice_surface_temp', temp)

    nc.define_2d_field('land_ice_area_fraction_retreat',
                       attrs={'units': '',
                              'long_name': 'mask where -ocean_kill cuts off ice'})
    nc.write('land_ice_area_fraction_retreat', okill)

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
    parser.description = "Creates a MISMIP boostrapping file for use with PISM."
    parser.add_option("-o", dest="output_filename",
                      help="output file")
    parser.add_option("-p", "--perturbation", dest="perturbation", type="string",
                      help="Perturbation at calving front (50m less at last cell (p1) or penultimate cell (p2))",
                      default="default")
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

    (opts, args) = parser.parse_args()


    if not opts.output_filename:
        output_filename = "ICESHELF_%d.nc" % opts.step
    else:
        output_filename = opts.output_filename

    print("Creating ICESHELF setup for step %s in %s..." % (
        opts.step, output_filename))

    config(opts.step)

    pism_bootstrap_file(output_filename,
                        opts.step,
                        opts.v0,
                        opts.H0,
                        calving_front=opts.calving_front,
                        N=opts.N,
                        p=opts.perturbation)

    print("done.")
