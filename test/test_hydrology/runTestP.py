#!/usr/bin/env python
from sys import argv, stderr, exit
import subprocess
import numpy as np
from argparse import ArgumentParser

try:
    from PISM import exactP
except:
    stderr.write("Please build PISM's Python bindings'.\n")
    exit(1)

try:
    from PISMNC import PISMDataset
except:
    subprocess.call("ln -sf ../../util/PISMNC.py", shell=True)
    from PISMNC import PISMDataset


def parse_options():
    stderr.write("reading options ...\n")

    parser = ArgumentParser()
    parser.description = "Test P (verification of '-hydrology distributed')."
    parser.add_argument("--pism_path", dest="PISM_PATH", default=".")
    parser.add_argument("--mpiexec", dest="MPIEXEC", default="")
    parser.add_argument("--Mx", dest="Mx",
                        help="Horizontal grid size. Default corresponds to a 1km grid.", type=int, default=51)
    parser.add_argument("--keep", dest="keep", action="store_true", help="Keep the generated PISM input file.")

    return parser.parse_args()


def generate_config():
    """Generates the config file with custom ice softness and hydraulic conductivity."""

    stderr.write("generating testPconfig.nc ...\n")

    nc = PISMDataset("testPconfig.nc", 'w')
    pism_overrides = nc.createVariable("pism_overrides", 'b')

    attrs = {
        "flow_law.isothermal_Glen.ice_softness": 3.1689e-24,
        "flow_law.isothermal_Glen.ice_softness_doc": "Pa-3 s-1; ice softness; NOT DEFAULT",

        "hydrology.hydraulic_conductivity": 1.0e-2 / (1000.0 * 9.81),
        "hydrology.hydraulic_conductivity_doc": "= k; NOT DEFAULT",

        "hydrology.regularizing_porosity": 0.01,
        "hydrology.regularizing_porosity_doc": "[pure]; phi_0 in notes",

        "hydrology.tillwat_max": 0.0,
        "hydrology.tillwat_max_doc": "m; turn off till water mechanism",

        "hydrology.thickness_power_in_flux": 1.0,
        "hydrology.thickness_power_in_flux_doc": "; = alpha in notes",

        "hydrology.gradient_power_in_flux": 2.0,
        "hydrology.gradient_power_in_flux_doc": "; = beta in notes",

        "hydrology.roughness_scale": 1.0,
        "hydrology.roughness_scale_doc": "m; W_r in notes; roughness scale",

        "basal_yield_stress.model": "constant",
        "basal_yield_stress.model_doc": "only the constant yield stress model works without till",

        "basal_yield_stress.constant.value": 1e6,
        "basal_yield_stress.constant.value_doc": "set default to 'high tauc'",
    }
    keys = list(attrs.keys())
    keys.sort()
    for k in keys:
        pism_overrides.setncattr(k, attrs[k])

    nc.close()


def report_drift(name, file1, file2, xx, yy, doshow=False):
    "Report on the difference between two files."
    nc1 = PISMDataset(file1)
    nc2 = PISMDataset(file2)

    var1 = nc1.variables[name]
    var2 = nc2.variables[name]
    diff = np.abs(np.squeeze(var1[:]) - np.squeeze(var2[:]))

    rr = np.sqrt(xx ** 2 + yy ** 2)
    diff[rr >= 0.89 * 25000.0] = 0.0

    if (doshow):
        import matplotlib.pyplot as plt
        plt.pcolormesh(xx, yy, diff)
        plt.axis('equal')
        plt.axis('tight')
        plt.colorbar()
        plt.show()

    #stderr.write("Drift in %s: average = %f, max = %f [%s]" % (name, np.average(diff), np.max(diff), var1.units) + "\n")
    return np.average(diff), np.max(diff)


def create_grid(Mx):
    Lx = 25.0e3  # outside L = 22.5 km

    x = np.linspace(-Lx, Lx, Mx)
    xx, yy = np.meshgrid(x, x)
    return x, x, xx, yy


def radially_outward(mag, x, y):
    """return components of a vector field  V(x,y)  which is radially-outward from
    the origin and has magnitude mag"""
    r = np.sqrt(x * x + y * y)
    if r == 0.0:
        return (0.0, 0.0)
    vx = mag * x / r
    vy = mag * y / r
    return (vx, vy)


def compute_sorted_radii(xx, yy):
    stderr.write("sorting radial variable ...\n")

    Mx = xx.shape[0]
    # create 1D array of tuples (r,j,k), sorted by r-value
    dtype = [('r', float), ('j', int), ('k', int)]
    rr = np.empty((Mx, Mx), dtype=dtype)

    for j in range(Mx):
        for k in range(Mx):
            rr[j, k] = (np.sqrt(xx[j, k] ** 2 + yy[j, k] ** 2), j, k)

    r = np.sort(rr.flatten(), order='r')

    return np.flipud(r)


def generate_pism_input(x, y, xx, yy):
    stderr.write("calling exactP() ...\n")

    EPS_ABS = 1.0e-12
    EPS_REL = 1.0e-15

    p = exactP(r[:]['r'], EPS_ABS, EPS_REL, 1)
    h_r, magvb_r, W_r, P_r = p.h, p.magvb, p.W, p.P

    stderr.write("creating gridded variables ...\n")
    # put on grid
    h = np.zeros_like(xx)
    W = np.zeros_like(xx)
    P = np.zeros_like(xx)

    magvb = np.zeros_like(xx)
    ussa = np.zeros_like(xx)
    vssa = np.zeros_like(xx)

    for n, pt in enumerate(r):
        j = pt['j']
        k = pt['k']
        h[j, k] = h_r[n]         # ice thickness in m
        magvb[j, k] = magvb_r[n]  # sliding speed in m s-1
        ussa[j, k], vssa[j, k] = radially_outward(magvb[j, k], xx[j, k], yy[j, k])
        W[j, k] = W_r[n]         # water thickness in m
        P[j, k] = P_r[n]         # water pressure in Pa

    stderr.write("creating inputforP.nc ...\n")

    nc = PISMDataset("inputforP.nc", 'w')
    nc.create_dimensions(x, y, time_dependent=True, use_time_bounds=True)

    nc.define_2d_field("thk", time_dependent=False,
                       attrs={"long_name": "ice thickness",
                              "units": "m",
                              "valid_min": 0.0,
                              "standard_name": "land_ice_thickness"})
    nc.define_2d_field("topg", time_dependent=False,
                       attrs={"long_name": "bedrock topography",
                              "units": "m",
                              "standard_name": "bedrock_altitude"})
    nc.define_2d_field("climatic_mass_balance", time_dependent=False,
                       attrs={"long_name": "climatic mass balance for -surface given",
                              "units": "kg m-2 year-1",
                              "standard_name": "land_ice_surface_specific_mass_balance"})
    nc.define_2d_field("ice_surface_temp", time_dependent=False,
                       attrs={"long_name": "ice surface temp (K) for -surface given",
                              "units": "Kelvin",
                              "valid_min": 0.0})
    nc.define_2d_field("bmelt", time_dependent=False,
                       attrs={"long_name": "basal melt rate",
                              "units": "m year-1",
                              "standard_name": "land_ice_basal_melt_rate"})

    nc.define_2d_field("bwat", time_dependent=False,
                       attrs={"long_name": "thickness of basal water layer",
                              "units": "m",
                              "valid_min": 0.0})
    nc.define_2d_field("bwp", time_dependent=False,
                       attrs={"long_name": "water pressure in basal water layer",
                              "units": "Pa",
                              "valid_min": 0.0})

    nc.define_2d_field("bc_mask", time_dependent=False,
                       attrs={"long_name": "if =1, apply u_ssa_bc and v_ssa_bc as sliding velocity"})
    nc.define_2d_field("u_ssa_bc", time_dependent=False,
                       attrs={"long_name": "x-component of prescribed sliding velocity",
                              "units": "m s-1"})
    nc.define_2d_field("v_ssa_bc", time_dependent=False,
                       attrs={"long_name": "y-component of prescribed sliding velocity",
                              "units": "m s-1"})

    Phi0 = 0.20                           # 20 cm/year basal melt rate
    T_surface = 260                       # ice surface temperature, K

    variables = {"topg": np.zeros_like(xx),
                 "climatic_mass_balance": np.zeros_like(xx),
                 "ice_surface_temp": np.ones_like(xx) + T_surface,
                 "bmelt": np.zeros_like(xx) + Phi0,
                 "thk": h,
                 "bwat": W,
                 "bwp": P,
                 "bc_mask": np.ones_like(xx),
                 "u_ssa_bc": ussa,
                 "v_ssa_bc": vssa}

    for name in list(variables.keys()):
        nc.write(name, variables[name])

    nc.history = subprocess.list2cmdline(argv)
    nc.close()
    stderr.write("NetCDF file %s written\n" % "inputforP.nc")


def run_pism(opts):
    stderr.write("Testing: Test P verification of '-hydrology distributed'.\n")

    cmd = "%s %s/pismr -config_override testPconfig.nc -i inputforP.nc -bootstrap -Mx %d -My %d -Mz 11 -Lz 4000 -hydrology distributed -report_mass_accounting -y 0.08333333333333 -max_dt 0.01 -no_mass -energy none -stress_balance ssa+sia -ssa_dirichlet_bc -o end.nc" % (
        opts.MPIEXEC, opts.PISM_PATH, opts.Mx, opts.Mx)

    stderr.write(cmd + "\n")
    subprocess.call(cmd, shell=True)

# high-res and parallel example:
#    ./runTestP.py --pism_path=../../build --mpiexec="mpiexec -n 4" --Mx=201
# example which should suffice for regression:
#    ./runTestP.py --pism_path=../../build --Mx=21


if __name__ == "__main__":

    opts = parse_options()

    x, y, xx, yy = create_grid(opts.Mx)

    r = compute_sorted_radii(xx, yy)

    generate_config()

    generate_pism_input(x, y, xx, yy)

    run_pism(opts)

    (bwatav, bwatmax) = report_drift("bwat", "inputforP.nc", "end.nc", xx, yy, doshow=False)
    (bwpav,  bwpmax) = report_drift("bwp",  "inputforP.nc", "end.nc", xx, yy, doshow=False)

    print("NUMERICAL ERRORS:")
    print("%d  %f  %f  %f  %f\n" % (opts.Mx, bwatav, bwatmax, bwpav, bwpmax))

    # cleanup:
    if opts.keep == False:
        subprocess.call("rm testPconfig.nc inputforP.nc end.nc", shell=True)
