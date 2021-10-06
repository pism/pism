"""This file contains very superficial tests of the PISM Python
wrappers. The goal is to be able to detect changes in the API
(function signatures, etc), not to test correctness.

Use with nose (https://pypi.python.org/pypi/nose/) and coverage.py
(https://pypi.python.org/pypi/coverage)

Run this to get a coverage report:

nosetests --with-coverage --cover-branches --cover-html --cover-package=PISM test/miscellaneous.py
"""

import PISM
import PISM.testing
import sys
import os
import numpy as np
from unittest import TestCase, SkipTest
from PISM.testing import filename

ctx = PISM.Context()
ctx.log.set_threshold(0)

def create_dummy_grid():
    "Create a dummy grid"
    ctx = PISM.Context()
    params = PISM.GridParameters(ctx.config)
    params.ownership_ranges_from_options(ctx.size)
    return PISM.IceGrid(ctx.ctx, params)


def context_test():
    "Test creating a new PISM context"
    ctx = PISM.Context()
    config = ctx.config
    us = ctx.unit_system
    EC = ctx.enthalpy_converter


def context_missing_attribute_test():
    "Test the handling of missing attributes"
    ctx = PISM.Context()
    try:
        config = ctx.foo        # there is no "foo", this should fail
        return False
    except AttributeError:
        return True


def create_grid_test():
    "Test the creation of the IceGrid object"
    grid1 = create_dummy_grid()

    grid2 = PISM.model.initGrid(PISM.Context(), 100e3, 100e3, 4000, 11, 11, 21, PISM.CELL_CORNER)


def algorithm_failure_exception_test():
    "Test the AlgorithmFailureException class"
    try:
        raise PISM.AlgorithmFailureException("no good reason")
        return False            # should not be reached
    except PISM.AlgorithmFailureException as e:
        print("calling e.reason(): ", e.reason())
        print("{}".format(e))
        return True


def printing_test():
    "Test verbPrintf"
    ctx = PISM.Context()
    PISM.verbPrintf(1, ctx.com, "hello %s!\n", "world")


def random_vec_test():
    "Test methods creating random fields"
    grid = PISM.IceGrid_Shallow(PISM.Context().ctx, 1e6, 1e6, 0, 0, 61, 31,
                                PISM.NOT_PERIODIC, PISM.CELL_CENTER)

    vec_scalar = PISM.vec.randVectorS(grid, 1.0)
    vec_vector = PISM.vec.randVectorV(grid, 2.0)

    vec_scalar_ghosted = PISM.vec.randVectorS(grid, 1.0, 2)
    vec_vector_ghosted = PISM.vec.randVectorV(grid, 2.0, 2)


def vec_metadata_test():
    "Test accessing IceModelVec metadata"
    grid = create_dummy_grid()

    vec_scalar = PISM.vec.randVectorS(grid, 1.0)

    m = vec_scalar.metadata()

    m.set_string("units", "kg")

    print(m.get_string("units"))


def vars_ownership_test():
    "Test passing IceModelVec ownership from Python to C++ (i.e. PISM)."
    grid = create_dummy_grid()
    variables = PISM.Vars()

    print("Adding 'thk'...")
    variables.add(PISM.model.createIceThicknessVec(grid))
    print("Returned from add_thk()...")

    print("Getting 'thk' from variables...")
    thk = variables.get("thk")
    print(thk)
    thk.begin_access()
    print("thickness at 0,0 is", thk[0, 0])
    thk.end_access()


def vec_access_test():
    "Test the PISM.vec.Access class and IceGrid::points, points_with_ghosts, coords"
    grid = create_dummy_grid()

    vec_scalar = PISM.vec.randVectorS(grid, 1.0)
    vec_scalar_ghosted = PISM.vec.randVectorS(grid, 1.0, 2)

    with PISM.vec.Access(comm=[vec_scalar_ghosted], nocomm=vec_scalar):
        for (i, j) in grid.points_with_ghosts():
            pass

    with PISM.vec.Access(comm=vec_scalar_ghosted, nocomm=[vec_scalar]):
        for (i, j) in grid.points():
            # do something
            pass

        for (i, j, x, y) in grid.coords():
            # do something with coordinates
            pass

    # try with nocomm=None
    with PISM.vec.Access(comm=vec_scalar_ghosted):
        pass


def create_modeldata_test():
    "Test creating the ModelData class"
    grid = create_dummy_grid()
    md = PISM.model.ModelData(grid)

    md2 = PISM.model.ModelData(grid, config=grid.ctx().config())


def grid_from_file_test():
    "Intiialize a grid from a file (test 1)"
    grid = create_dummy_grid()

    enthalpy = PISM.model.createEnthalpyVec(grid)
    enthalpy.set(80e3)

    file_name = filename("grid_from_file")
    try:
        pio = PISM.util.prepare_output(file_name)

        enthalpy.write(pio)

        pio = PISM.File(grid.com, file_name, PISM.PISM_NETCDF3, PISM.PISM_READONLY)

        grid2 = PISM.IceGrid.FromFile(ctx.ctx, pio, "enthalpy", PISM.CELL_CORNER)
    finally:
        os.remove(file_name)

def grid_from_file_test_2():
    "Intiialize a grid from a file (test 2)"
    grid = create_dummy_grid()

    enthalpy = PISM.model.createEnthalpyVec(grid)
    enthalpy.set(80e3)

    file_name = filename("grid_from_file")
    try:
        pio = PISM.util.prepare_output(file_name)

        enthalpy.write(pio)

        grid2 = PISM.IceGrid.FromFile(ctx.ctx, file_name, ["enthalpy"], PISM.CELL_CORNER)
    finally:
        os.remove(file_name)

def create_special_vecs_test():
    "Test helpers used to create standard PISM fields"
    grid = create_dummy_grid()

    usurf = PISM.model.createIceSurfaceVec(grid)

    thk = PISM.model.createIceThicknessVec(grid)

    sea_level = PISM.model.createSeaLevelVec(grid)

    usurfstore = PISM.model.createIceSurfaceStoreVec(grid)

    thkstore = PISM.model.createIceThicknessStoreVec(grid)

    bed = PISM.model.createBedrockElevationVec(grid)

    tauc = PISM.model.createYieldStressVec(grid)

    strainheat = PISM.model.createStrainHeatingVec(grid)

    u, v, w = PISM.model.create3DVelocityVecs(grid)

    hardav = PISM.model.createAveragedHardnessVec(grid)

    enthalpy = PISM.model.createEnthalpyVec(grid)

    age = PISM.model.createAgeVec(grid)

    bmr = PISM.model.createBasalMeltRateVec(grid)

    tillphi = PISM.model.createTillPhiVec(grid)

    basal_water = PISM.model.createBasalWaterVec(grid)

    gl_mask = PISM.model.createGroundingLineMask(grid)

    vel = PISM.model.create2dVelocityVec(grid)

    taudx = PISM.model.createDrivingStressXVec(grid)

    taudy = PISM.model.createDrivingStressYVec(grid)

    vel_misfit_weight = PISM.model.createVelocityMisfitWeightVec(grid)

    cbar = PISM.model.createCBarVec(grid)

    mask = PISM.model.createIceMaskVec(grid)

    bcmask = PISM.model.createBCMaskVec(grid)

    no_model_mask = PISM.model.createNoModelMaskVec(grid)

    zeta_fixed_mask = PISM.model.createZetaFixedMaskVec(grid)

    lon = PISM.model.createLongitudeVec(grid)

    lat = PISM.model.createLatitudeVec(grid)

    # test ModelVecs.add()
    modeldata = PISM.model.ModelData(grid)
    vecs = modeldata.vecs

    vecs.add(mask)

    print(vecs)
    # test getattr
    vecs.mask

    return True


def pism_vars_test():
    """Test adding fields to and getting them from pism::Vars."""
    grid = create_dummy_grid()

    v = grid.variables()

    v.add(PISM.model.createIceThicknessVec(grid))

    # test getting by short name
    print(v.get("thk").metadata().get_string("units"))

    # test getting by standard name
    print(v.get("land_ice_thickness").metadata().get_string("units"))


def modelvecs_test():
    "Test the ModelVecs class"

    grid = create_dummy_grid()

    mask = PISM.model.createIceMaskVec(grid)
    mask.set(PISM.MASK_GROUNDED)

    modeldata = PISM.model.ModelData(grid)
    vecs = modeldata.vecs

    vecs.add(mask, "ice_mask", writing=True)

    # use the default name, no writing
    vecs.add(PISM.model.createIceThicknessVec(grid))

    try:
        vecs.add(mask, "ice_mask")
        return False
    except RuntimeError:
        # should fail: mask was added already
        pass

    # get a field:
    print("get() method: ice mask: ", vecs.get("ice_mask").metadata().get_string("long_name"))

    print("dot notation: ice mask: ", vecs.ice_mask.metadata().get_string("long_name"))

    try:
        vecs.invalid
        return False
    except AttributeError:
        # should fail
        pass

    try:
        vecs.get("invalid")
        return False
    except RuntimeError:
        # should fail
        pass

    # test __repr__
    print(vecs)

    # test has()
    print("Has thickness?", vecs.has("thickness"))

    # test markForWriting
    vecs.markForWriting("ice_mask")

    vecs.markForWriting(mask)

    vecs.markForWriting("thk")

    # test write()
    output_file = filename("test_ModelVecs")
    try:
        pio = PISM.util.prepare_output(output_file)
        pio.close()

        vecs.write(output_file)

        # test writeall()
        vecs.writeall(output_file)
    finally:
        os.remove(output_file)

def sia_test():
    "Test the PISM.sia module"
    ctx = PISM.Context()
    params = PISM.GridParameters(ctx.config)
    params.Lx = 1e5
    params.Ly = 1e5
    params.Lz = 1000
    params.Mx = 100
    params.My = 100
    params.Mz = 11
    params.registration = PISM.CELL_CORNER
    params.periodicity = PISM.NOT_PERIODIC
    params.ownership_ranges_from_options(ctx.size)
    grid = PISM.IceGrid(ctx.ctx, params)

    enthalpyconverter = PISM.EnthalpyConverter(ctx.config)

    mask = PISM.model.createIceMaskVec(grid)
    mask.set(PISM.MASK_GROUNDED)

    thk = PISM.model.createIceThicknessVec(grid)
    thk.set(1000.0)

    surface = PISM.model.createIceSurfaceVec(grid)
    surface.set(1000.0)

    bed = PISM.model.createBedrockElevationVec(grid)
    bed.set(0.0)

    enthalpy = PISM.model.createEnthalpyVec(grid)
    enthalpy.set(enthalpyconverter.enthalpy(270.0, 0.0, 0.0))

    modeldata = PISM.model.ModelData(grid)
    modeldata.setPhysics(enthalpyconverter)

    vecs = grid.variables()

    fields = [thk, surface, mask, bed, enthalpy]

    for field in fields:
        vecs.add(field)

    vel_sia = PISM.sia.computeSIASurfaceVelocities(modeldata)


def util_test():
    "Test the PISM.util module"
    grid = create_dummy_grid()

    output_file = filename("test_pism_util")
    try:
        pio = PISM.File(grid.com, output_file, PISM.PISM_NETCDF3, PISM.PISM_READWRITE_MOVE)
        pio.close()

        PISM.util.writeProvenance(output_file)
        PISM.util.writeProvenance(output_file, message="history string")

        PISM.util.fileHasVariable(output_file, "data")
    finally:
        os.remove(output_file)

    # Test PISM.util.Bunch
    b = PISM.util.Bunch(a=1, b="string")
    b.update(c=3.0)

    print(b.a, b["b"], "b" in b, b)


def logging_test():
    "Test the PISM.logging module"
    grid = create_dummy_grid()

    import PISM.logging as L

    log_filename = filename("log")
    try:
        PISM.File(grid.com, log_filename, PISM.PISM_NETCDF3, PISM.PISM_READWRITE_MOVE)
        c = L.CaptureLogger(log_filename)

        L.clear_loggers()

        L.add_logger(L.print_logger)
        L.add_logger(c)

        L.log("log message\n", L.kError)

        L.logError("error message\n")

        L.logWarning("warning message\n")

        L.logMessage("log message (again)\n")

        L.logDebug("debug message\n")

        L.logPrattle("prattle message\n")

        c.write()                   # default arguments
        c.readOldLog()
    finally:
        os.remove(log_filename)

    log_filename = filename("other_log")
    try:
        PISM.File(grid.com, log_filename, PISM.PISM_NETCDF3, PISM.PISM_READWRITE_MOVE)
        c.write(log_filename, "other_log")  # non-default arguments
    finally:
        os.remove(log_filename)


def column_interpolation_test(plot=False):
    """Test ColumnInterpolation by interpolating from the coarse grid to the
    fine grid and back."""
    import numpy as np

    Lz = 1000.0
    Mz = 41

    def z_quadratic(Mz, Lz):
        "Compute levels of a quadratic coarse grid."
        result = np.zeros(Mz)
        z_lambda = 4.0
        for k in range(Mz - 1):
            zeta = float(k) / (Mz - 1)
            result[k] = Lz * ((zeta / z_lambda) * (1.0 + (z_lambda - 1.0) * zeta))
        result[Mz - 1] = Lz
        return result

    def fine_grid(z_coarse):
        "Compute levels of the fine grid corresponding to a given coarse grid."
        Lz = z_coarse[-1]
        dz = np.min(np.diff(z_coarse))
        Mz = int(np.ceil(Lz / dz) + 1)
        dz = Lz / (Mz - 1.0)
        result = np.zeros(Mz)
        for k in range(1, Mz):
            result[k] = z_coarse[0] + k * dz

        return result

    def test_quadratic_interp():
        z_coarse = z_quadratic(Mz, Lz)
        f_coarse = (z_coarse / Lz) ** 2
        z_fine = fine_grid(z_coarse)

        print("Testing quadratic interpolation")
        return test_interp(z_coarse, f_coarse, z_fine, "Quadratic interpolation")

    def test_linear_interp():
        z_coarse = np.linspace(0, Lz, Mz)
        f_coarse = (z_coarse / Lz) ** 2
        z_fine = fine_grid(z_coarse)

        print("Testing linear interpolation")
        return test_interp(z_coarse, f_coarse, z_fine, "Linear interpolation")

    def test_interp(z, f, z_fine, title):
        interp = PISM.ColumnInterpolation(z, z_fine)

        f_fine = interp.coarse_to_fine(f, interp.Mz_fine())

        f_fine_numpy = np.interp(z_fine, z, f)

        f_roundtrip = interp.fine_to_coarse(f_fine)

        if plot:
            plt.figure()
            plt.plot(z, f, 'o-', label="original coarse-grid data")
            plt.plot(z_fine, f_fine, 'o-', label="interpolated onto the fine grid")
            plt.plot(z, f_roundtrip, 'o-', label="interpolated back onto the coarse grid")
            plt.plot(z, f_roundtrip - f, 'o-', label="difference after the roundtrip")
            plt.legend(loc="best")
            plt.title(title)
            plt.grid(True)


        delta = np.linalg.norm(f - f_roundtrip, ord=1)
        delta_numpy = np.linalg.norm(f_fine - f_fine_numpy, ord=1)
        print("norm1(fine_to_coarse(coarse_to_fine(f)) - f) = %f" % delta)
        print("norm1(PISM - NumPy) = %f" % delta_numpy)

        return delta, delta_numpy

    if plot:
        import pylab as plt

    linear_delta, linear_delta_numpy = test_linear_interp()

    quadratic_delta, _ = test_quadratic_interp()

    if plot:
        plt.show()

    if (linear_delta > 1e-12 or linear_delta_numpy > 1e-12 or quadratic_delta > 1e-3):
        return False
    return True


def pism_join_test():
    "Test PISM.join()"
    assert PISM.join(["one", "two"], ':') == "one:two"


def pism_split_test():
    "Test PISM.split()"
    assert PISM.split("one,two,three", ',') == ("one", "two", "three")


def pism_ends_with_test():
    "Test PISM.ends_with()"
    assert PISM.ends_with("foo.nc", ".nc") == True
    assert PISM.ends_with("foo.nc and more text", ".nc") == False
    assert PISM.ends_with("short_string", "longer_suffix") == False


def linear_interpolation_test(plot=False):
    "Test linear interpolation code used to regrid fields"
    import numpy as np

    M_in = 11
    M_out = 101
    a = 0.0
    b = 10.0
    padding = 1.0
    x_input = np.linspace(a, b, M_in)
    x_output = np.sort(((b + padding) - (a - padding)) * np.random.rand(M_out) + (a - padding))

    def F(x):
        return x * 2.0 + 5.0

    values = F(x_input)

    i = PISM.Interpolation(PISM.LINEAR, x_input, x_output)

    F_interpolated = i.interpolate(values)

    F_desired = F(x_output)
    F_desired[x_output < a] = F(a)
    F_desired[x_output > b] = F(b)

    if plot:
        import pylab as plt

        plt.plot(x_output, F_interpolated, 'o-', color='blue', label="interpolated result")
        plt.plot(x_output, F_desired, 'x-', color='green', label="desired result")
        plt.plot(x_input, values, 'o-', color='red', label="input")
        plt.grid(True)
        plt.legend(loc="best")
        plt.show()

    assert np.max(np.fabs(F_desired - F_interpolated)) < 1e-16


def pism_context_test():
    "Test creating and using a C++-level Context"

    com = PISM.PETSc.COMM_WORLD
    system = PISM.UnitSystem("")

    logger = PISM.Logger(com, 2)

    config = PISM.DefaultConfig(com, "pism_config", "-config", system)
    config.init_with_default(logger)

    EC = PISM.EnthalpyConverter(config)

    time = PISM.Time(com, config, logger, system)

    ctx = PISM.cpp.Context(com, system, config, EC, time, logger, "greenland")

    print(ctx.com().Get_size())
    print(ctx.config().get_number("constants.standard_gravity"))
    print(ctx.enthalpy_converter().L(273.15))
    print(ctx.time().current())
    print(PISM.convert(ctx.unit_system(), 1, "km", "m"))
    print(ctx.prefix())


def check_flow_law(factory, flow_law_name, EC, stored_data):
    factory.set_default(flow_law_name)
    law = factory.create()

    depth = 2000
    gs = 1e-3
    sigma = [1e4, 5e4, 1e5, 1.5e5]

    T_pa = [-30, -5, 0, 0]
    omega = [0.0, 0.0, 0.0, 0.005]

    assert len(T_pa) == len(omega)

    p = EC.pressure(depth)
    Tm = EC.melting_temperature(p)

    data = []
    print("  Flow table for %s" % law.name())
    print("| Sigma        | Temperature  | Omega        | Flow factor  |")
    print("|--------------+--------------+--------------+--------------|")
    for S in sigma:
        for Tpa, O in zip(T_pa, omega):
            T = Tm + Tpa
            E = EC.enthalpy(T, O, p)
            F = law.flow(S, E, p, gs)
            data.append(F)

            print("| %e | %e | %e | %e |" % (S, T, O, F))
    print("|--------------+--------------+--------------+--------------|")
    print("")

    data = np.array(data)

    assert np.max(np.fabs(data - stored_data)) < 1e-16


def flowlaw_test():
    data = {}
    data["arr"] = [3.91729503e-18, 6.42803396e-17, 1.05746828e-16, 1.05746828e-16,
                   9.79323757e-17, 1.60700849e-15, 2.64367070e-15, 2.64367070e-15,
                   3.91729503e-16, 6.42803396e-15, 1.05746828e-14, 1.05746828e-14,
                   8.81391381e-16, 1.44630764e-14, 2.37930363e-14, 2.37930363e-14]
    data["arrwarm"] = [1.59798478e-19, 1.04360343e-16, 3.30653997e-16, 3.30653997e-16,
                       3.99496194e-18, 2.60900856e-15, 8.26634991e-15, 8.26634991e-15,
                       1.59798478e-17, 1.04360343e-14, 3.30653997e-14, 3.30653997e-14,
                       3.59546574e-17, 2.34810771e-14, 7.43971492e-14, 7.43971492e-14]
    data["gk"] = [1.1636334595808724e-16, 6.217445758362754e-15, 2.5309103327753672e-14,
                  2.5309103327753672e-14, 2.5947947614616463e-16, 2.0065832524499375e-14,
                  9.158056141786197e-14, 9.158056141786197e-14, 4.493111202368685e-16,
                  3.469816186746473e-14, 1.6171243121742907e-13, 1.6171243121742907e-13,
                  7.12096200221403e-16, 4.879162291119208e-14, 2.2895389865988545e-13, 2.2895389865988545e-13]
    data["gpbld"] = [4.65791754e-18, 1.45114704e-16, 4.54299921e-16, 8.66009225e-16,
                     1.16447938e-16, 3.62786761e-15, 1.13574980e-14, 2.16502306e-14,
                     4.65791754e-16, 1.45114704e-14, 4.54299921e-14, 8.66009225e-14,
                     1.04803145e-15, 3.26508084e-14, 1.02217482e-13, 1.94852076e-13]
    data["hooke"] = [5.26775897e-18, 2.12325906e-16, 5.32397091e-15, 5.32397091e-15,
                     1.31693974e-16, 5.30814764e-15, 1.33099273e-13, 1.33099273e-13,
                     5.26775897e-16, 2.12325906e-14, 5.32397091e-13, 5.32397091e-13,
                     1.18524577e-15, 4.77733287e-14, 1.19789346e-12, 1.19789346e-12]
    data["isothermal_glen"] = [3.16890000e-16, 3.16890000e-16, 3.16890000e-16, 3.16890000e-16,
                               7.92225000e-15, 7.92225000e-15, 7.92225000e-15, 7.92225000e-15,
                               3.16890000e-14, 3.16890000e-14, 3.16890000e-14, 3.16890000e-14,
                               7.13002500e-14, 7.13002500e-14, 7.13002500e-14, 7.13002500e-14]
    data["pb"] = [4.65791754e-18, 1.45114704e-16, 4.54299921e-16, 4.54299921e-16,
                  1.16447938e-16, 3.62786761e-15, 1.13574980e-14, 1.13574980e-14,
                  4.65791754e-16, 1.45114704e-14, 4.54299921e-14, 4.54299921e-14,
                  1.04803145e-15, 3.26508084e-14, 1.02217482e-13, 1.02217482e-13]

    ctx = PISM.context_from_options(PISM.PETSc.COMM_WORLD, "flowlaw_test")
    EC = ctx.enthalpy_converter()
    factory = PISM.FlowLawFactory("stress_balance.sia.", ctx.config(), EC)

    for flow_law_name, data in data.items():
        check_flow_law(factory, flow_law_name, EC, np.array(data))


def ssa_trivial_test():
    "Test the SSA solver using a trivial setup."

    context = PISM.Context()

    L = 50.e3  # // 50km half-width
    H0 = 500  # // m
    dhdx = 0.005  # // pure number, slope of surface & bed
    nu0 = PISM.util.convert(30.0, "MPa year", "Pa s")
    tauc0 = 1.e4  # // 1kPa

    class TrivialSSARun(PISM.ssa.SSAExactTestCase):
        def _initGrid(self):
            self.grid = PISM.IceGrid.Shallow(context.ctx, L, L, 0, 0,
                                             self.Mx, self.My, PISM.CELL_CORNER, PISM.NOT_PERIODIC)

        def _initPhysics(self):
            self.modeldata.setPhysics(context.enthalpy_converter)

        def _initSSACoefficients(self):
            self._allocStdSSACoefficients()
            self._allocateBCs()

            vecs = self.modeldata.vecs

            vecs.land_ice_thickness.set(H0)
            vecs.surface_altitude.set(H0)
            vecs.bedrock_altitude.set(0.0)
            vecs.tauc.set(tauc0)

            # zero Dirichler B.C. everywhere
            vecs.vel_bc.set(0.0)
            vecs.vel_bc_mask.set(1.0)

        def _initSSA(self):
            # The following ensure that the strength extension is used everywhere
            se = self.ssa.strength_extension
            se.set_notional_strength(nu0 * H0)
            se.set_min_thickness(4000 * 10)

            # For the benefit of SSAFD on a non-periodic grid
            self.config.set_flag("ssa.compute_surface_gradient_inward", True)

        def exactSolution(self, i, j, x, y):
            return [0, 0]

    output_file = filename("ssa_trivial")
    try:
        Mx = 11
        My = 11
        test_case = TrivialSSARun(Mx, My)
        test_case.run(output_file)
    finally:
        os.remove(output_file)

def epsg_test():
    "Test EPSG to CF conversion."
    l = PISM.StringLogger(PISM.PETSc.COMM_WORLD, 2)
    system = PISM.Context().unit_system

    # test supported formats
    for template in ["{epsg}:{code}",
                     "+init={epsg}:{code}",
                     "+units=m +init={epsg}:{code}",
                     "+init={epsg}:{code} +units=m"]:
        for epsg in ["EPSG", "epsg"]:
            for code in [3413, 3031, 3057, 5936, 26710]:
                string = template.format(epsg=epsg, code=code)
                print("Trying {}".format(string))
                l.reset()
                v = PISM.epsg_to_cf(system, string)
                v.report_to_stdout(l, 2)
                print(l.get())
                print("done.")

    # test that unsupported codes trigger an exception
    try:
        v = PISM.epsg_to_cf(system, "+init=epsg:3032")
        raise AssertionError("should fail with 3032: only 3413 and 3031 are supported")
    except RuntimeError as e:
        print("unsupported codes trigger exceptions: {}".format(e))

    # test that an invalid PROJ string (e.g. an EPSG code is not a
    # number) triggers an exception
    try:
        v = PISM.epsg_to_cf(system, "+init=epsg:not-a-number +units=m")
        raise AssertionError("an invalid PROJ string failed to trigger an exception")
    except RuntimeError as e:
        print("invalid codes trigger exceptions: {}".format(e))

def regridding_test():
    "Test 2D regridding: same input and target grids."
    import numpy as np

    ctx = PISM.Context()
    params = PISM.GridParameters(ctx.config)
    params.Mx = 3
    params.My = 3
    params.ownership_ranges_from_options(1)

    grid = PISM.IceGrid(ctx.ctx, params)

    thk1 = PISM.model.createIceThicknessVec(grid)
    thk2 = PISM.model.createIceThicknessVec(grid)
    x = grid.x()
    x_min = np.min(x)
    x_max = np.max(x)
    y = grid.y()
    y_min = np.min(y)
    y_max = np.max(y)
    with PISM.vec.Access(nocomm=[thk1]):
        for (i, j) in grid.points():
            F_x = (x[i] - x_min) / (x_max - x_min)
            F_y = (y[j] - y_min) / (y_max - y_min)
            thk1[i, j] = (F_x + F_y) / 2.0

    file_name = filename("thickness")
    try:
        thk1.dump(file_name)

        thk2.regrid(file_name, critical=True)

        with PISM.vec.Access(nocomm=[thk1, thk2]):
            for (i, j) in grid.points():
                v1 = thk1[i, j]
                v2 = thk2[i, j]
                if np.abs(v1 - v2) > 1e-12:
                    raise AssertionError("mismatch at {},{}: {} != {}".format(i, j, v1, v2))

    finally:
        os.remove(file_name)

def interpolation_weights_test():
    "Test 2D interpolation weights."

    def interp2d(grid, F, x, y):
        i_left, i_right, j_bottom, j_top = grid.compute_point_neighbors(x, y)
        w = grid.compute_interp_weights(x, y)

        i = [i_left, i_right, i_right, i_left]
        j = [j_bottom, j_bottom, j_top, j_top]

        result = 0.0
        for k in range(4):
            result += w[k] * F[j[k], i[k]]

        return result

    Mx = 100
    My = 200
    Lx = 20
    Ly = 10

    grid = PISM.IceGrid_Shallow(PISM.Context().ctx,
                                Lx, Ly, 0, 0, Mx, My,
                                PISM.CELL_CORNER,
                                PISM.NOT_PERIODIC)

    x = grid.x()
    y = grid.y()
    X, Y = np.meshgrid(x, y)
    Z = 2 * X + 3 * Y

    N = 1000
    np.random.seed(1)
    x_pts = np.random.rand(N) * (2 * Lx) - Lx
    y_pts = np.random.rand(N) * (2 * Ly) - Ly
    # a linear function should be recovered perfectly
    exact = 2 * x_pts + 3 * y_pts

    result = np.array([interp2d(grid, Z, x_pts[k], y_pts[k]) for k in range(N)])

    np.testing.assert_almost_equal(result, exact)


def vertical_extrapolation_during_regridding_test():
    "Test extrapolation in the vertical direction"
    # create a grid with 11 levels, 1000m thick
    ctx = PISM.Context()
    params = PISM.GridParameters(ctx.config)
    params.Lx = 1e5
    params.Ly = 1e5
    params.Mx = 3
    params.My = 3
    params.Mz = 11
    params.Lz = 1000
    params.registration = PISM.CELL_CORNER
    params.periodicity = PISM.NOT_PERIODIC
    params.ownership_ranges_from_options(ctx.size)

    z = np.linspace(0, params.Lz, params.Mz)
    params.z[:] = z

    grid = PISM.IceGrid(ctx.ctx, params)

    # create an IceModelVec that uses this grid
    v = PISM.IceModelVec3(grid, "test", PISM.WITHOUT_GHOSTS, grid.z())
    v.set(0.0)

    # set a column
    with PISM.vec.Access(nocomm=[v]):
        v.set_column(1, 1, z)

    # save to a file
    file_name = filename("regridding")
    try:
        v.dump(file_name)

        # create a taller grid (to 2000m):
        params.Lz = 2000
        params.Mz = 41
        z_tall = np.linspace(0, params.Lz, params.Mz)
        params.z[:] = z_tall

        tall_grid = PISM.IceGrid(ctx.ctx, params)

        # create an IceModelVec that uses this grid
        v_tall = PISM.IceModelVec3(tall_grid, "test", PISM.WITHOUT_GHOSTS, tall_grid.z())

        # Try regridding without extrapolation. This should fail.
        try:
            ctx.ctx.log().disable()
            v_tall.regrid(file_name, PISM.CRITICAL)
            ctx.ctx.log().enable()
            raise AssertionError("Should not be able to regrid without extrapolation")
        except RuntimeError as e:
            pass

        # allow extrapolation during regridding
        ctx.config.set_flag("grid.allow_extrapolation", True)

        # regrid from test.nc
        ctx.ctx.log().disable()
        v_tall.regrid(file_name, PISM.CRITICAL)
        ctx.ctx.log().enable()

        # get a column
        with PISM.vec.Access(nocomm=[v_tall]):
            column = np.array(v_tall.get_column(1, 1))

        # compute the desired result
        desired = np.r_[np.linspace(0, 1000, 21), np.zeros(20) + 1000]

        # compare
        np.testing.assert_almost_equal(column, desired)
    finally:
        os.remove(file_name)

class PrincipalStrainRates(TestCase):
    def u_exact(self, x, y):
        "Velocity field for testing"
        return (np.cos(y) + np.sin(x),
                np.sin(y) + np.cos(x))

    def eps_exact(self, x, y):
        "Principal strain rates corresponding to u_exact."
        u_x = np.cos(x)
        u_y = - np.sin(y)

        v_x = - np.sin(x)
        v_y = np.cos(y)

        A   = 0.5 * (u_x + v_y)
        B   = 0.5 * (u_x - v_y)
        Dxy = 0.5 * (v_x + u_y)
        q   = np.sqrt(B**2 + Dxy**2);

        return (A + q, A - q)

    def create_grid(self, Mx):
        My = Mx + 10            # a non-square grid
        return PISM.IceGrid.Shallow(self.ctx.ctx,
                                    2*np.pi, 2*np.pi,
                                    0, 0, int(Mx), int(My), PISM.CELL_CENTER, PISM.NOT_PERIODIC)

    def create_velocity(self, grid):
        velocity = PISM.IceModelVec2V(grid, "bar", PISM.WITH_GHOSTS)
        with PISM.vec.Access(nocomm=velocity):
            for (i, j) in grid.points():
                u, v = self.u_exact(grid.x(i), grid.y(j))
                velocity[i, j].u = u
                velocity[i, j].v = v
        velocity.update_ghosts()

        return velocity

    def create_cell_type(self, grid):
        cell_type = PISM.IceModelVec2CellType(grid, "cell_type", PISM.WITH_GHOSTS)
        cell_type.set(PISM.MASK_GROUNDED)
        cell_type.update_ghosts()

        return cell_type

    def error(self, Mx):
        grid = self.create_grid(Mx)

        velocity = self.create_velocity(grid)
        cell_type = self.create_cell_type(grid)
        strain_rates = PISM.IceModelVec3(grid, "strain_rates",
                                         PISM.WITHOUT_GHOSTS, 2)

        PISM.compute_2D_principal_strain_rates(velocity, cell_type, strain_rates)
        rates = strain_rates.numpy()

        e1 = rates[:,:,0]
        e2 = rates[:,:,1]

        # compute e1, e2 using formulas
        xx, yy = np.meshgrid(grid.x(), grid.y())
        e1_exact, e2_exact = self.eps_exact(xx, yy)

        delta_e1 = np.max(np.fabs(e1 - e1_exact))
        delta_e2 = np.max(np.fabs(e2 - e2_exact))

        return Mx, delta_e1, delta_e2

    def setUp(self):
        self.ctx = PISM.Context()

    def test_principal_strain_rates(self):
        "Test principal strain rate computation"
        errors = np.array([self.error(M) for M in [10, 20, 40]])

        Ms = errors[:,0]
        p_e1 = np.polyfit(np.log(1.0 / Ms), np.log(errors[:,1]), 1)
        p_e2 = np.polyfit(np.log(1.0 / Ms), np.log(errors[:,2]), 1)

        assert p_e1[0] > 1.8
        assert p_e2[0] > 1.8

    def tearDown(self):
        pass

def test_trapezoid_integral():
    "Linear integration weights"
    x = [0.0, 0.5, 1.0, 2.0]
    y = [2.0, 2.5, 3.0, 2.0]

    def f(a, b):
        weights = PISM.integration_weights(x, PISM.LINEAR, a, b)

        result = 0.0
        for k, w in weights.items():
            result += w * y[k]
        return result

    assert f(0, 2) == 5.0
    assert f(0, 0.5) + f(0.5, 1) + f(1, 1.5) + f(1.5, 2) == f(0, 2)
    assert f(0, 0.5) + f(0.5, 1.5) + f(1.5, 2) == f(0, 2)
    assert f(0, 0.25) + f(0.25, 1.75) + f(1.75, 2) == f(0, 2)

    # constant extrapolation:
    assert f(-2, -1) == 1 * y[0]
    assert f(-2, 0) == 2 * y[0]
    assert f(2, 5) == 3 * y[-1]
    assert f(3, 4) == 1 * y[-1]


def test_piece_wise_constant_integral():
    "Linear integration weights"
    x = [0.0, 0.5, 1.0, 2.0]
    y = [2.0, 2.5, 3.0, 2.0]

    def f(a, b):
        weights = PISM.integration_weights(x, PISM.PIECEWISE_CONSTANT, a, b)

        result = 0.0
        for k, w in weights.items():
            result += w * y[k]
        return result

    np.testing.assert_almost_equal(f(0, 2), 5.25)
    np.testing.assert_almost_equal(f(0.25, 0.35), 0.2)
    np.testing.assert_almost_equal(f(0, 0.5) + f(0.5, 1) + f(1, 1.5) + f(1.5, 2), f(0, 2))
    np.testing.assert_almost_equal(f(0, 0.5) + f(0.5, 1.5) + f(1.5, 2), f(0, 2))
    np.testing.assert_almost_equal(f(0, 0.25) + f(0.25, 1.75) + f(1.75, 2), f(0, 2))

    # constant extrapolation:
    np.testing.assert_almost_equal(f(-2, -1), 1 * y[0])
    np.testing.assert_almost_equal(f(-2, 0), 2 * y[0])
    np.testing.assert_almost_equal(f(-2, 0.5), 2 * y[0] + 0.5 * 2)
    np.testing.assert_almost_equal(f(2, 5), 3 * y[-1])
    np.testing.assert_almost_equal(f(3, 4), 1 * y[-1])

def test_interpolation_other():
    x = [0.0, 10.0]

    try:
        PISM.Interpolation(PISM.PIECEWISE_CONSTANT, [2, 1], [2, 1])
        assert False, "failed to detect non-increasing times"
    except RuntimeError as e:
        print(e)
        pass

    try:
        W = PISM.integration_weights(x, PISM.NEAREST, 1, 2)
        assert False, "failed to detect an unsupported interpolation type"
    except RuntimeError as e:
        print(e)
        pass

    try:
        W = PISM.integration_weights(x, PISM.LINEAR, 1, 1)
        assert False, "failed to catch a >= b"
    except RuntimeError as e:
        print(e)
        pass

    I = PISM.Interpolation(PISM.LINEAR, [0, 1, 2], [0.5, 1.5])
    np.testing.assert_almost_equal(I.left(), [0, 1])
    np.testing.assert_almost_equal(I.right(), [1, 2])
    np.testing.assert_almost_equal(I.alpha(), [0.5, 0.5])

def test_nearest_neighbor():
    "Nearest neighbor interpolation"
    x = [-1, 1]
    y = [0, 1]

    xx = np.linspace(-0.9, 0.9, 10)
    yy = np.ones_like(xx) * (xx > 0)

    zz = PISM.Interpolation(PISM.NEAREST, x, xx).interpolate(y)

    np.testing.assert_almost_equal(yy, zz)

class AgeModel(TestCase):
    def setUp(self):
        self.output_file = filename("age")
        self.grid = self.create_dummy_grid()

    def create_dummy_grid(self):
        "Create a dummy grid"
        params = PISM.GridParameters(ctx.config)
        params.ownership_ranges_from_options(ctx.size)
        return PISM.IceGrid(ctx.ctx, params)

    def test_age_model_runs(self):
        "Check if AgeModel runs"
        ice_thickness = PISM.model.createIceThicknessVec(self.grid)

        u = PISM.IceModelVec3(self.grid, "u", PISM.WITHOUT_GHOSTS, self.grid.z())
        v = PISM.IceModelVec3(self.grid, "v", PISM.WITHOUT_GHOSTS, self.grid.z())
        w = PISM.IceModelVec3(self.grid, "w", PISM.WITHOUT_GHOSTS, self.grid.z())

        ice_thickness.set(4000.0)
        u.set(0.0)
        v.set(0.0)
        w.set(0.0)

        model = PISM.AgeModel(self.grid, None)
        input_options = PISM.process_input_options(ctx.com, ctx.config)
        model.init(input_options)

        inputs = PISM.AgeModelInputs(ice_thickness, u, v, w)

        dt = PISM.util.convert(1, "years", "seconds")

        model.update(0, dt, inputs)

        model.age().dump(self.output_file)

    def tearDown(self):
        os.remove(self.output_file)

def checksum_test():
    "Check if a small change in an IceModelVec affects checksum() output"
    grid = PISM.testing.shallow_grid(Mx=101, My=201)

    v = PISM.IceModelVec2S(grid, "dummy", PISM.WITHOUT_GHOSTS)
    v.set(1e15)

    old_checksum = v.checksum()

    with PISM.vec.Access(nocomm=v):
        for (i, j) in grid.points():
            if i == 0 and j == 0:
                v[i, j] += 1

    assert old_checksum != v.checksum()

class ForcingOptions(TestCase):
    def setUp(self):
        # store current configuration parameters
        self.config = PISM.DefaultConfig(ctx.com, "pism_config", "-config", ctx.unit_system)
        self.config.init_with_default(ctx.log)
        self.config.import_from(ctx.config)

        self.filename = filename("forcing-options-input")
        PISM.util.prepare_output(self.filename)

    def test_without_file(self):
        "ForcingOptions: xxx.file is not set"
        ctx.config.set_string("input.file", self.filename)

        opt = PISM.ForcingOptions(ctx.ctx, "surface.given")

        assert opt.filename == self.filename
        assert opt.periodic == False

    def test_with_file(self):
        "ForcingOptions: xxx.file is set"
        ctx.config.set_string("surface.given.file", self.filename)

        opt = PISM.ForcingOptions(ctx.ctx, "surface.given")

        assert opt.filename == self.filename
        assert opt.periodic == False

    def test_without_file_and_without_input_file(self):
        "ForcingOptions: xxx.file is not set and -i is not set"
        try:
            opt = PISM.ForcingOptions(ctx.ctx, "surface.given")
            assert False, "failed to stop with an error message"
        except RuntimeError:
            pass

    def tearDown(self):
        # reset configuration parameters
        ctx.config.import_from(self.config)

        os.remove(self.filename)

def test_bq2():
    "2-point boundary quadrature"
    raise SkipTest("not implemented")
    dx = 1.0
    dy = dx
    Q = PISM.BoundaryQuadrature2(1, 1, 1)
    for q in [0, 1]:
        for s in [0, 1, 2, 3]:
            psi = [Q.germ(s, q, k).val for k in [0, 1]]
            assert sum(psi) == 1.0, "side = {}, q = {}, psi = {}".format(s, q, psi)

class ScalarForcing(TestCase):
    def setUp(self):
        suffix = filename("-scalar-forcing-")
        self.filename              = "input" + suffix
        self.filename_decreasing   = "decreasing" + suffix
        self.filename_wrong_bounds = "wrong-bounds" + suffix
        self.filename_2d           = "2d" + suffix
        self.filename_1            = "one-value" + suffix
        self.filename_no_bounds    = "no-bounds" + suffix

        def f(x):
            return np.sin(2 * np.pi / 12 * x)

        dt = 1.0
        self.times = (0.5 + np.arange(12)) * dt
        self.values = f(self.times)

        self.times_long = (0.5 + np.arange(24)) * dt
        self.values_long = f(self.times_long)

        self.ts = np.linspace(0 + dt, 24 - dt, 1001)

        time_bounds = np.vstack((self.times - 0.5 * dt, self.times + 0.5 * dt)).T.flatten()

        ctx.time.set_start(0)
        ctx.time.set_end(4)

        # good input
        PISM.testing.create_scalar_forcing(self.filename,
                                           "delta_T", "Kelvin",
                                           self.values, self.times, time_bounds)

        # non-increasing times
        PISM.testing.create_scalar_forcing(self.filename_decreasing,
                                           "delta_T", "Kelvin",
                                           [0, 1], [1, 0], time_bounds=[1, 2, 0, 1])

        # wrong time bounds
        PISM.testing.create_scalar_forcing(self.filename_wrong_bounds,
                                           "delta_T", "Kelvin",
                                           [0, 1], [0, 1], time_bounds=[0, 2, 2, 3])
        # invalid number of dimensions
        grid = create_dummy_grid()
        PISM.testing.create_forcing(grid,
                                    self.filename_2d,
                                    "delta_T",
                                    "Kelvin",
                                    [0, 1],
                                    "seconds since 1-1-1",
                                    times=[0, 1])

        # only one value
        PISM.testing.create_scalar_forcing(self.filename_1,
                                           "delta_T", "Kelvin",
                                           [1], [0], time_bounds=[0, 4])

        # no time bounds
        PISM.testing.create_scalar_forcing(self.filename_no_bounds,
                                           "delta_T", "Kelvin",
                                           [1], [0], time_bounds=None)

    def tearDown(self):
        try:
            os.remove(self.filename)
            os.remove(self.filename_decreasing)
            os.remove(self.filename_wrong_bounds)
            os.remove(self.filename_2d)
            os.remove(self.filename_1)
            os.remove(self.filename_no_bounds)
        except:
            pass

    def create(self, filename=None, periodic=False):
        if filename is not None:
            return PISM.ScalarForcing(ctx.ctx,
                                      filename,
                                      "delta_T",
                                      "K",
                                      "K",
                                      "temperature offsets", periodic)

        return PISM.ScalarForcing(ctx.ctx,
                                  "surface.delta_T",
                                  "delta_T",
                                  "K",
                                  "K",
                                  "temperature offsets")

    def test_file_not_set(self):
        try:
            ctx.config.set_string("surface.delta_T.file", "")
            self.create()
            assert False, "failed to stop if prefix.file is empty"
        except RuntimeError:
            pass

    def test_decreasing_time(self):
        try:
            self.create(self.filename_decreasing, False)
            assert False, "failed to stop if times are non-increasing"
        except RuntimeError as e:
            print(e)
            pass

    def test_wrong_bounds(self):
        try:
            self.create(self.filename_wrong_bounds, False)
            assert False, "failed to stop if time bounds are wrong"
        except RuntimeError as e:
            print(e)
            pass

    def test_invalid_dimensions(self):
        try:
            self.create(self.filename_2d, False)
            assert False, "failed to stop if the variable is not scalar"
        except RuntimeError as e:
            print(e)
            pass

    def test_no_time_bounds(self):
        try:
            self.create(self.filename_no_bounds, False)
            assert False, "failed to stop if time bounds are missing"
        except RuntimeError:
            pass

    def test_run_duration(self):
        end = ctx.time.end()
        try:
            ctx.time.set_end(100)
            self.create(self.filename, False)
            assert False, "failed to stop because forcing does not span model time"
        except RuntimeError as e:
            print(e)
            pass
        finally:
            ctx.time.set_end(end)

    def test_periodic_interpolation(self):
        F = self.create(self.filename, True)

        Y_numpy = np.interp(self.ts, self.times_long, self.values_long)
        Y = [F.value(t) for t in self.ts]

        np.testing.assert_almost_equal(Y, Y_numpy)

    def test_interpolation(self):
        ctx.config.set_string("surface.delta_T.file", self.filename)
        ctx.config.set_flag("surface.delta_T.periodic", False)

        F = self.create()

        # two outside and one inside
        T = [-1, 13, 5.123]

        Y_numpy = np.interp(T, self.times, self.values)
        Y = [F.value(t) for t in T]

        np.testing.assert_almost_equal(Y, Y_numpy)

    def test_one_value(self):
        F = self.create(self.filename_1, False)

        T = [1, 2]

        Y_ref = [1, 1]
        Y = [F.value(t) for t in T]

        np.testing.assert_almost_equal(Y, Y_ref)

    def test_average_scalar_forcing(self):
        """Test ScalarForcing::average(t, dt)"""

        def compute_average(times, time_bounds, values, t, dt, periodic=False):
            filename = PISM.testing.filename("forcing-")
            try:
                PISM.testing.create_scalar_forcing(filename, "delta_T", "Kelvin",
                                                   values, times, time_bounds)

                F = self.create(filename, periodic)

                return F.average(t, dt)
            finally:
                os.remove(filename)

        times       = [0.5, 1.5, 2.5, 3.5]
        time_bounds = [0, 1, 1, 2, 2, 3, 3, 4]
        values      = [1, 1, -1, -1]

        def f(a, b):
            return compute_average(times, time_bounds, values,
                                   a, b - a, periodic=True)

        def g(a, b):
            return compute_average(times, time_bounds, values,
                                   a, b - a, periodic=False)

        # both points are outside on the left
        np.testing.assert_equal(g(-2, -1), 1.0)
        # both are outside on the right
        np.testing.assert_equal(g(4.5, 5.5), -1.0)
        # one too far left, one too far right
        np.testing.assert_equal(g(-0.5, 4.5), 0.0)
        # same sub-interval
        np.testing.assert_equal(g(0, 0.25), 1.0)
        # dt == 0.0
        np.testing.assert_equal(g(0, 0), 1.0)

        try:
            np.testing.assert_equal(g(1, 0), 1.0)
            assert False, "failed to catch negative dt"
        except RuntimeError:
            pass

        # one full period
        np.testing.assert_equal(f(0, 4), 0.0)
        # two periods
        np.testing.assert_equal(f(0, 8), 0.0)
        # one period, but starting from 4
        np.testing.assert_equal(f(4, 8), 0.0)
        # part of a period
        np.testing.assert_equal(f(4, 6), 0.75)
        # part of a period (check if shifting affects results)
        np.testing.assert_equal(f(4, 6), f(0, 2))

def thickness_calving_test():
    "Test the time-dependent thickness calving threshold"

    grid = PISM.testing.shallow_grid(Mx=5, My=5)

    geometry = PISM.Geometry(grid)

    # set ice thickness
    H = geometry.ice_thickness
    with PISM.vec.Access(H):
        for i,j in grid.points():
            if i >= 3 and j >= 3:
                H[i, j] = 0
            elif i > 2 and j == 2:
                H[i, j] = 100
            elif i >= 1 and j >= 1:
                H[i, j] = 150
            else:
                H[i, j] = 200

    # set bed elevation low enough to ensure that all ice is floating
    geometry.bed_elevation.set(-1000)

    geometry.ensure_consistency(0.0)

    config = PISM.Context().config

    filename = PISM.testing.filename("threshold_thickness_")
    day = 86400.0
    times = [0, 30, 60]
    time_bounds = [-15, 15, 45, 75]
    time_units = "days since 1-1-1"
    values = [100, 150, 200]
    try:
        # create the time-dependent thickness threshold
        PISM.testing.create_forcing(grid,
                                    filename,
                                    "thickness_calving_threshold",
                                    "m",
                                    values,
                                    time_units,
                                    times=times,
                                    time_bounds=time_bounds)

        config.set_string("calving.thickness_calving.file", filename)
        calving = PISM.CalvingAtThickness(grid)
        calving.init()

        # run the calving model
        calving.update(0*day, 30*day, geometry.cell_type, geometry.ice_thickness)

        ref1 = np.array([[200., 200., 200., 200., 200.],
                         [200., 150., 150., 150., 150.],
                         [200., 150., 150.,   0.,   0.],
                         [200., 150., 150.,   0.,   0.],
                         [200., 150., 150.,   0.,   0.]])

        np.testing.assert_almost_equal(geometry.ice_thickness.numpy(), ref1)

        calving.update(30*day, 60*day, geometry.cell_type, geometry.ice_thickness)

        ref2 = np.array([[200., 200., 200., 200., 200.],
                         [200., 150., 150.,   0.,   0.],
                         [200., 150.,   0.,   0.,   0.],
                         [200., 150.,   0.,   0.,   0.],
                         [200., 150.,   0.,   0.,   0.]])

        np.testing.assert_almost_equal(geometry.ice_thickness.numpy(), ref2)

    finally:
        # clean up
        os.remove(filename)

        # restore default values
        config.set_string("calving.thickness_calving.file", "")
