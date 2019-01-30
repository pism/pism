"""This file contains very superficial tests of the PISM Python
wrappers. The goal is to be able to detect changes in the API
(function signatures, etc), not to test correctness.

Use with nose (https://pypi.python.org/pypi/nose/) and coverage.py
(https://pypi.python.org/pypi/coverage)

Run this to get a coverage report:

nosetests --with-coverage --cover-branches --cover-html --cover-package=PISM test/nosetests.py
"""

import PISM
import sys
import numpy as np


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
    "Intiialize a grid from a file"
    grid = create_dummy_grid()

    enthalpy = PISM.model.createEnthalpyVec(grid)
    enthalpy.set(80e3)

    output_file = "test_grid_from_file.nc"
    pio = PISM.util.prepare_output(output_file)

    enthalpy.write(pio)

    pio = PISM.PIO(grid.com, "netcdf3", output_file, PISM.PISM_READONLY)

    grid2 = PISM.IceGrid.FromFile(grid.ctx(), pio, "enthalpy", PISM.CELL_CORNER)


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
    output_file = "test_ModelVecs.nc"
    pio = PISM.util.prepare_output(output_file)
    pio.close()

    vecs.write(output_file)

    # test writeall()
    vecs.writeall(output_file)


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

    output_file = "test_pism_util.nc"
    pio = PISM.PIO(grid.com, "netcdf3", output_file, PISM.PISM_READWRITE_MOVE)
    pio.close()

    PISM.util.writeProvenance(output_file)
    PISM.util.writeProvenance(output_file, message="history string")

    PISM.util.fileHasVariable(output_file, "data")

    # Test PISM.util.Bunch
    b = PISM.util.Bunch(a=1, b="string")
    b.update(c=3.0)

    print(b.a, b["b"], "b" in b, b)


def logging_test():
    "Test the PISM.logging module"
    grid = create_dummy_grid()

    import PISM.logging as L

    PISM.PIO(grid.com, "netcdf3", "log.nc", PISM.PISM_READWRITE_MOVE)
    c = L.CaptureLogger("log.nc")

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

    PISM.PIO(grid.com, "netcdf3", "other_log.nc", PISM.PISM_READWRITE_MOVE)
    c.write("other_log.nc", "other_log")  # non-default arguments


def column_interpolation_test(plot=False):
    """Test ColumnInterpolation by interpolating from the coarse grid to the
    fine grid and back."""
    import numpy as np
    import pylab as plt

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

        def plot():
            plt.figure()
            plt.plot(z, f, 'o-', label="original coarse-grid data")
            plt.plot(z_fine, f_fine, 'o-', label="interpolated onto the fine grid")
            plt.plot(z, f_roundtrip, 'o-', label="interpolated back onto the coarse grid")
            plt.plot(z, f_roundtrip - f, 'o-', label="difference after the roundtrip")
            plt.legend(loc="best")
            plt.title(title)
            plt.grid(True)

        if plot:
            plot()

        delta = np.linalg.norm(f - f_roundtrip, ord=1)
        delta_numpy = np.linalg.norm(f_fine - f_fine_numpy, ord=1)
        print("norm1(fine_to_coarse(coarse_to_fine(f)) - f) = %f" % delta)
        print("norm1(PISM - NumPy) = %f" % delta_numpy)

        return delta, delta_numpy

    linear_delta, linear_delta_numpy = test_linear_interp()

    quadratic_delta, _ = test_quadratic_interp()

    if plot:
        plt.show()

    if (linear_delta > 1e-12 or
            linear_delta_numpy > 1e-12 or
            quadratic_delta > 1e-3):
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

    i = PISM.LinearInterpolation(x_input, x_output)

    F_interpolated = i.interpolate(values)

    F_desired = F(x_output)
    F_desired[x_output < a] = F(a)
    F_desired[x_output > b] = F(b)

    if plot:
        import pylab as plt

        plt.hold(True)
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

    time = PISM.Time(config, "360_day", system)

    ctx = PISM.cpp.Context(com, system, config, EC, time, logger, "greenland")

    print(ctx.com().Get_size())
    print(ctx.config().get_double("constants.standard_gravity"))
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
            vecs.bc_mask.set(1.0)

        def _initSSA(self):
            # The following ensure that the strength extension is used everywhere
            se = self.ssa.strength_extension
            se.set_notional_strength(nu0 * H0)
            se.set_min_thickness(4000 * 10)

            # For the benefit of SSAFD on a non-periodic grid
            self.config.set_boolean("ssa.compute_surface_gradient_inward", True)

        def exactSolution(self, i, j, x, y):
            return [0, 0]

    Mx = 11
    My = 11
    test_case = TrivialSSARun(Mx, My)
    test_case.run("ssa_trivial.nc")


def epsg_test():
    "Test EPSG to CF conversion."
    l = PISM.StringLogger(PISM.PETSc.COMM_WORLD, 2)
    system = PISM.Context().unit_system

    # test supported EPSG codes
    for code in [3413, 3031]:
        print("Trying code {}".format(code))
        l.reset()
        # +init at the beginning
        v = PISM.epsg_to_cf(system, "+init=epsg:%d" % code)
        # +init not at the beginning of the string
        v = PISM.epsg_to_cf(system, "+units=m +init=epsg:%d" % code)
        # +init followed by more options
        v = PISM.epsg_to_cf(system, "+init=epsg:%d +units=m" % code)
        v.report_to_stdout(l, 2)
        print(l.get())
        print("done.")

    # test that unsupported codes trigger an exception
    try:
        v = PISM.epsg_to_cf(system, "+init=epsg:3032")
        raise AssertionError("should fail with 3032: only 3413 and 3031 are supported")
    except RuntimeError as e:
        print("unsupported codes trigger exceptions: {}".format(e))

    # test that an invalid PROJ.4 string (e.g. an EPSG code is not a
    # number) triggers an exception
    try:
        v = PISM.epsg_to_cf(system, "+init=epsg:not-a-number +units=m")
        # raise AssertionError("an invalid PROJ.4 string failed to trigger an exception")
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
    thk1.dump("thk1.nc")

    thk2.regrid("thk1.nc", critical=True)

    with PISM.vec.Access(nocomm=[thk1, thk2]):
        for (i, j) in grid.points():
            v1 = thk1[i, j]
            v2 = thk2[i, j]
            if np.abs(v1 - v2) > 1e-12:
                raise AssertionError("mismatch at {},{}: {} != {}".format(i, j, v1, v2))

    import os
    os.remove("thk1.nc")


def netcdf_string_attribute_test():
    "Test reading a NetCDF-4 string attribute."
    import os

    basename = "string_attribute_test"
    attribute = "string attribute"

    def setup():
        cdl = """
netcdf string_attribute_test {
  string :string_attribute = "%s" ;
  :text_attribute = "%s" ;
}
""" % (attribute, attribute)
        with open(basename + ".cdl", "w") as f:
            f.write(cdl)

        os.system("ncgen -4 %s.cdl" % basename)

    def teardown():
        # remove the temporary file
        os.remove(basename + ".nc")
        os.remove(basename + ".cdl")

    def compare(backend):
        try:
            pio = PISM.PIO(PISM.PETSc.COMM_WORLD, backend, basename + ".nc", PISM.PISM_READONLY)
        except:
            # Don't fail if backend creation failed: PISM may not have
            # been compiled with parallel I/O enabled.
            return

        read_string = pio.get_att_text("PISM_GLOBAL", "string_attribute")
        read_text = pio.get_att_text("PISM_GLOBAL", "text_attribute")

        # check that written and read strings are the same
        print("written string: '%s'" % attribute)
        print("read string:    '%s'" % read_string)
        print("read text:      '%s'" % read_text)
        assert read_string == attribute
        assert read_text == attribute

    setup()

    compare("netcdf3")
    compare("netcdf4_parallel")

    teardown()


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
    v = PISM.IceModelVec3()
    v.create(grid, "test", PISM.WITHOUT_GHOSTS)
    v.set(0.0)

    # set a column
    with PISM.vec.Access(nocomm=[v]):
        v.set_column(1, 1, z)

    # save to a file
    v.dump("test.nc")

    # create a taller grid (to 2000m):
    params.Lz = 2000
    params.Mz = 41
    z_tall = np.linspace(0, params.Lz, params.Mz)
    params.z[:] = z_tall

    tall_grid = PISM.IceGrid(ctx.ctx, params)

    # create an IceModelVec that uses this grid
    v_tall = PISM.IceModelVec3()
    v_tall.create(tall_grid, "test", PISM.WITHOUT_GHOSTS)

    # Try regridding without extrapolation. This should fail.
    try:
        ctx.ctx.log().disable()
        v_tall.regrid("test.nc", PISM.CRITICAL)
        ctx.ctx.log().enable()
        raise AssertionError("Should not be able to regrid without extrapolation")
    except RuntimeError as e:
        pass

    # allow extrapolation during regridding
    ctx.config.set_boolean("grid.allow_extrapolation", True)

    # regrid from test.nc
    ctx.ctx.log().disable()
    v_tall.regrid("test.nc", PISM.CRITICAL)
    ctx.ctx.log().enable()

    # get a column
    with PISM.vec.Access(nocomm=[v_tall]):
        column = np.array(v_tall.get_column_vector(1, 1))

    # compute the desired result
    desired = np.r_[np.linspace(0, 1000, 21), np.zeros(20) + 1000]

    # compare
    np.testing.assert_almost_equal(column, desired)

    # clean up
    import os
    os.remove("test.nc")
