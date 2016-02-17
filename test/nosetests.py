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

    grid2 = PISM.model.initGrid(PISM.Context(), 100e3, 100e3, 4000, 11, 11, 21, PISM.NOT_PERIODIC)


def algorithm_failure_exception_test():
    "Test the AlgorithmFailureException class"
    try:
        raise PISM.AlgorithmFailureException("no good reason")
        return False            # should not be reached
    except PISM.AlgorithmFailureException as e:
        print "calling e.reason(): ", e.reason()
        print "{}".format(e)
        return True


def printing_test():
    "Test verbPrintf"
    ctx = PISM.Context()
    PISM.verbPrintf(1, ctx.com, "hello %s!\n", "world")


def random_vec_test():
    "Test methods creating random fields"
    grid = create_dummy_grid()

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

    print m.get_string("units")


def vars_ownership_test():
    "Test passing IceModelVec ownership from Python to C++ (i.e. PISM)."
    grid = create_dummy_grid()
    variables = PISM.Vars()

    print "Adding 'thk'..."
    variables.add(PISM.model.createIceThicknessVec(grid))
    print "Returned from add_thk()..."

    print "Getting 'thk' from variables..."
    thk = variables.get("thk")
    print thk
    thk.begin_access()
    print "thickness at 0,0 is", thk[0, 0]
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


def toproczero_test():
    "Test communication to processor 0"
    grid = create_dummy_grid()

    vec_scalar = PISM.vec.randVectorS(grid, 1.0)
    vec_vector = PISM.vec.randVectorV(grid, 2.0)

    tz = PISM.vec.ToProcZero(grid)
    array_scalar_0 = tz.communicate(vec_scalar)

    tz2 = PISM.vec.ToProcZero(grid, dof=2, dim=2)
    array_vector_0 = tz2.communicate(vec_vector)

    try:
        tz3 = PISM.vec.ToProcZero(grid, dof=2, dim=3)
        return False
    except NotImplementedError:
        # 3D fields are not supported (yet)
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
    pio = PISM.PIO(grid.com, "netcdf3")
    pio.open(output_file, PISM.PISM_READWRITE_MOVE)
    PISM.define_time(pio, grid.ctx().config().get_string("time_dimension_name"),
                     grid.ctx().config().get_string("calendar"),
                     grid.ctx().time().units_string(),
                     grid.ctx().unit_system())
    PISM.append_time(pio,
                     grid.ctx().config().get_string("time_dimension_name"),
                     grid.ctx().time().current())
    pio.close()

    enthalpy.write(output_file)

    grid2 = PISM.IceGrid.FromFile(grid.ctx(), output_file, "enthalpy", PISM.NOT_PERIODIC)

def create_special_vecs_test():
    "Test helpers used to create standard PISM fields"
    grid = create_dummy_grid()

    usurf = PISM.model.createIceSurfaceVec(grid)

    thk = PISM.model.createIceThicknessVec(grid)

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

    cell_area = PISM.model.createCellAreaVec(grid)

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

    print vecs
    # test getattr
    vecs.mask

    return True


def options_test():
    "Test command-line option handling"
    ctx = PISM.Context()

    o = PISM.PETSc.Options()

    M = PISM.optionsInt("-M", "description", default=100)
    M = PISM.optionsInt("-M", "description", default=None)

    S = PISM.optionsString("-S", "description", default="string")
    S = PISM.optionsString("-S", "description", default=None)

    R = PISM.optionsReal("-R", "description", default=1.5)
    R = PISM.optionsReal("-R", "description", default=None)

    o.setValue("-B", "on")
    B = PISM.optionsFlag("-B", "description", default=False)
    B = PISM.optionsFlag("B", "description", default=False)
    B = PISM.optionsFlag("-B", "description", default=None)

    o.setValue("-no_C", "on")
    C = PISM.optionsFlag("C", "description", default=None)

    D = PISM.optionsFlag("D", "description", default=None)
    D = PISM.optionsFlag("D", "description", default=True)

    o.setValue("-no_D", "on")
    o.setValue("-D", "on")
    try:
        # should throw RuntimeError
        D = PISM.optionsFlag("D", "description", default=None)
        return False
    except RuntimeError:
        pass

    o.setValue("-IA", "1,2,3")
    IA = PISM.optionsIntArray("-IA", "description", default=[1, 2])
    IA = PISM.optionsIntArray("-IA", "description", default=None)
    IA2 = PISM.optionsIntArray("-IA2", "description", default=None)
    IA2 = PISM.optionsIntArray("-IA2", "description", default=[1, 2])

    o.setValue("-RA", "1,2,3")
    RA = PISM.optionsRealArray("-RA", "description", default=[2, 3])
    RA = PISM.optionsRealArray("-RA", "description", default=None)
    RA2 = PISM.optionsRealArray("-RA2", "description", default=[2, 3])
    RA2 = PISM.optionsRealArray("-RA2", "description", default=None)

    o.setValue("-SA", "1,2,3")
    SA = PISM.optionsStringArray("-SA", "description", default="one,two")
    SA = PISM.optionsStringArray("-SA", "description", default=None)
    SA2 = PISM.optionsStringArray("-SA2", "description", default="two,three")
    SA2 = PISM.optionsStringArray("-SA2", "description", default=None)

    M = PISM.optionsList("-L", "description", choices="one,two", default="one")
    M = PISM.optionsList("-L", "description", choices="one,two", default=None)


def pism_vars_test():
    """Test adding fields to and getting them from pism::Vars."""
    grid = create_dummy_grid()

    v = grid.variables()

    v.add(PISM.model.createIceThicknessVec(grid))

    # test getting by short name
    print v.get("thk").metadata().get_string("units")

    # test getting by standard name
    print v.get("land_ice_thickness").metadata().get_string("units")


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
    print "get() method: ice mask: ", vecs.get("ice_mask").metadata().get_string("long_name")

    print "dot notation: ice mask: ", vecs.ice_mask.metadata().get_string("long_name")

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
    print vecs

    # test has()
    print "Has thickness?", vecs.has("thickness")

    # test markForWriting
    vecs.markForWriting("ice_mask")

    vecs.markForWriting(mask)

    vecs.markForWriting("thk")

    # test write()
    output_file = "test_ModelVecs.nc"
    pio = PISM.PIO(grid.com, "netcdf3")
    pio.open(output_file, PISM.PISM_READWRITE_MOVE)
    PISM.define_time(pio, grid.ctx().config().get_string("time_dimension_name"),
                     grid.ctx().config().get_string("calendar"),
                     grid.ctx().time().units_string(),
                     grid.ctx().unit_system())
    PISM.append_time(pio,
                     grid.ctx().config().get_string("time_dimension_name"),
                     grid.ctx().time().current())
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
    pio = PISM.PIO(grid.com, "netcdf3")
    pio.open(output_file, PISM.PISM_READWRITE_MOVE)
    pio.close()

    PISM.util.writeProvenance(output_file)
    PISM.util.writeProvenance(output_file, message="history string")

    PISM.util.fileHasVariable(output_file, "data")

    # Test PISM.util.Bunch
    b = PISM.util.Bunch(a=1, b="string")
    b.update(c=3.0)

    print b.a, b["b"], b.has_key("b"), b


def logging_test():
    "Test the PISM.logging module"
    grid = create_dummy_grid()
    pio = PISM.PIO(grid.com, "netcdf3")

    import PISM.logging as L

    pio.open("log.nc", PISM.PISM_READWRITE_MOVE)
    pio.close()
    c = L.CaptureLogger("log.nc")

    L.clear_loggers()

    L.add_logger(L.print_logger)
    L.add_logger(c)

    PISM.setVerbosityLevel(2)

    L.log("log message\n", L.kError)

    L.logError("error message\n")

    L.logWarning("warning message\n")

    L.logMessage("log message (again)\n")

    L.logDebug("debug message\n")

    L.logPrattle("prattle message\n")

    c.write()                   # default arguments
    c.readOldLog()

    pio.open("other_log.nc", PISM.PISM_READWRITE_MOVE)
    pio.close()
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
        for k in xrange(Mz - 1):
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

        print "Testing quadratic interpolation"
        return test_interp(z_coarse, f_coarse, z_fine, "Quadratic interpolation")

    def test_linear_interp():
        z_coarse = np.linspace(0, Lz, Mz)
        f_coarse = (z_coarse / Lz) ** 2
        z_fine = fine_grid(z_coarse)

        print "Testing linear interpolation"
        return test_interp(z_coarse, f_coarse, z_fine, "Linear interpolation")

    def test_interp(z, f, z_fine, title):
        interp = PISM.ColumnInterpolation(z, z_fine)

        f_fine = interp.coarse_to_fine(f, interp.Mz_fine())

        f_fine_numpy = np.interp(z_fine, z, f)

        f_roundtrip = interp.fine_to_coarse(f_fine)

        def plot():
            plt.figure()
            plt.hold(True)
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
        print "norm1(fine_to_coarse(coarse_to_fine(f)) - f) = %f" % delta
        print "norm1(PISM - NumPy) = %f" % delta_numpy

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

    print ctx.com().Get_size()
    print ctx.config().get_double("standard_gravity")
    print ctx.enthalpy_converter().L(273.15)
    print ctx.time().current()
    print PISM.convert(ctx.unit_system(), 1, "km", "m")
    print ctx.prefix()

def flowlaw_test():
    ctx = PISM.context_from_options(PISM.PETSc.COMM_WORLD, "flowlaw_test")
    EC = ctx.enthalpy_converter()
    ff = PISM.FlowLawFactory("sia_", ctx.config(), EC)
    law = ff.create()

    TpaC = [-30, -5, 0, 0]
    depth = 2000
    gs = 1e-3
    omega = [0.0, 0.0, 0.0, 0.005]
    sigma = [1e4, 5e4, 1e5, 1.5e5]

    p = EC.pressure(depth)
    Tm = EC.melting_temperature(p)

    print "flow law:   \"%s\"" % law.name()
    print "pressure = %9.3e Pa = (hydrostatic at depth %7.2f m)" % (p, depth)
    print "flowtable:"
    print "  (dev stress)   (abs temp) (liq frac) =   (flow)"

    for i in range(4):
        for j in range(4):
            T = Tm + TpaC[j]
            E = EC.enthalpy(T, omega[j], p)
            flowcoeff = law.flow(sigma[i], E, p, gs)
            print "    %10.2e   %10.3f  %9.3f = %10.6e" % (sigma[i], T, omega[j], flowcoeff)

def gpbld3_flow_test():
    "Test the optimized version of GPBLD."
    ctx = PISM.context_from_options(PISM.PETSc.COMM_WORLD, "GPBLD3_test")
    EC = ctx.enthalpy_converter()
    gpbld = PISM.GPBLD("sia_", ctx.config(), EC)
    gpbld3 = PISM.GPBLD3("sia_", ctx.config(), EC)

    import numpy as np
    N = 11
    TpaC = np.linspace(-30, 0, N)
    depth = np.linspace(0, 4000, N)
    omega = np.linspace(0, 0.02, N)
    sigma = [1e4, 5e4, 1e5, 1.5e5]

    gs = 1e-3

    for d in depth:
        p = EC.pressure(d)
        Tm = EC.melting_temperature(p)
        for Tpa in TpaC:
            T = Tm + Tpa
            for o in omega:
                if T >= Tm:
                    E = EC.enthalpy(T, o, p)
                else:
                    E = EC.enthalpy(T, 0.0, p)
                for s in sigma:
                    regular = gpbld.flow(s, E, p, gs)
                    optimized = gpbld3.flow(s, E, p, gs)
                    assert np.fabs(regular - optimized) / regular < 2e-14

def gpbld3_hardness_test():
    "Test the hardness implementation in the optimized version of GPBLD."
    ctx = PISM.context_from_options(PISM.PETSc.COMM_WORLD, "GPBLD3_test")
    EC = ctx.enthalpy_converter()
    gpbld = PISM.GPBLD("sia_", ctx.config(), EC)
    gpbld3 = PISM.GPBLD3("sia_", ctx.config(), EC)

    import numpy as np
    N = 11
    TpaC = np.linspace(-30, 0, N)
    depth = np.linspace(0, 4000, N)
    omega = np.linspace(0, 0.02, N)

    for d in depth:
        p = EC.pressure(d)
        Tm = EC.melting_temperature(p)
        for Tpa in TpaC:
            T = Tm + Tpa
            for o in omega:
                if T >= Tm:
                    E = EC.enthalpy(T, o, p)
                else:
                    E = EC.enthalpy(T, 0.0, p)

                regular = gpbld.hardness(E, p)
                optimized = gpbld3.hardness(E, p)
                assert np.fabs(regular - optimized) / regular < 4e-15


def gpbld3_error_report():
    """Print max. absolute and relative difference between GPBLD and
    GPBLD3. Uses 101*101*101*101 samples in a "reasonable" range of
    pressure-adjusted temperatures, depth, water fraction, and
    effective stress. This takes about 15 minutes to complete.
    """
    ctx = PISM.context_from_options(PISM.PETSc.COMM_WORLD, "GPBLD3_test")
    EC = ctx.enthalpy_converter()
    gpbld = PISM.GPBLD("sia_", ctx.config(), EC)
    gpbld3 = PISM.GPBLD3("sia_", ctx.config(), EC)

    import numpy as np
    N = 31
    TpaC = np.linspace(-30, 0, N)
    depth = np.linspace(0, 5000, N)
    omega = np.linspace(0, 0.02, N)
    sigma = np.linspace(0, 5e5, N)

    gs = 1e-3

    max_difference = 0.0
    max_rel_difference = 0.0

    for d in depth:
        p = EC.pressure(d)
        Tm = EC.melting_temperature(p)
        for Tpa in TpaC:
            T = Tm + Tpa
            for o in omega:
                if T >= Tm:
                    E = EC.enthalpy(T, o, p)
                else:
                    E = EC.enthalpy(T, 0.0, p)
                for s in sigma:
                    regular = gpbld.flow(s, E, p, gs)
                    optimized = gpbld3.flow(s, E, p, gs)
                    max_difference = max(np.fabs(regular - optimized), max_difference)
                    if regular > 0.0:
                        max_rel_difference = max(np.fabs(regular - optimized) / regular,
                                                 max_rel_difference)

    print "%d (%e) samples" % (N**4, N**4)
    print "max difference", max_difference
    print "max relative difference", max_rel_difference

def ssa_trivial_test():
    "Test the SSA solver using a trivial setup."

    context = PISM.Context()
    unit_system = context.unit_system

    L = 50.e3  # // 50km half-width
    H0 = 500  # // m
    dhdx = 0.005  # // pure number, slope of surface & bed
    nu0 = PISM.convert(unit_system, 30.0, "MPa year", "Pa s")
    tauc0 = 1.e4  # // 1kPa

    class TrivialSSARun(PISM.ssa.SSAExactTestCase):
        def _initGrid(self):
            self.grid = PISM.IceGrid.Shallow(PISM.Context().ctx, L, L, 0, 0,
                                             self.Mx, self.My, PISM.NONE)

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
            self.config.set_boolean("compute_surf_grad_inward_ssa", True)

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
        print "Trying code {}".format(code)
        l.reset()
        # +init at the beginning
        v = PISM.epsg_to_cf(system, "+init=epsg:%d" % code)
        # +init not at the beginning of the string
        v = PISM.epsg_to_cf(system, "+units=m +init=epsg:%d" % code)
        # +init followed by more options
        v = PISM.epsg_to_cf(system, "+init=epsg:%d +units=m" % code)
        v.report_to_stdout(l, 2)
        print l.get(),
        print "done."

    # test that unsupported codes trigger an exception
    try:
        v = PISM.epsg_to_cf(system, "+init=epsg:3032")
        raise AssertionError("should fail with 3032: only 3413 and 3031 are supported")
    except RuntimeError as e:
        print "unsupported codes trigger exceptions: {}".format(e)

    # test that an invalid PROJ.4 string (e.g. an EPSG code is not a
    # number) triggers an exception
    try:
        v = PISM.epsg_to_cf(system, "+init=epsg:not-a-number +units=m")
        # raise AssertionError("an invalid PROJ.4 string failed to trigger an exception")
    except RuntimeError as e:
        print "invalid codes trigger exceptions: {}".format(e)

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
            v1 = thk1[i,j]
            v2 = thk2[i,j]
            if np.abs(v1 - v2) > 1e-12:
                raise AssertionError("mismatch at {},{}: {} != {}".format(i, j, v1, v2))

    import os
    os.remove("thk1.nc")
def po_constant_test():
    """Test that the basal melt rate computed by ocean::Constant is the
    same regardless of whether it is set using
    ocean_sub_shelf_heat_flux_into_ice or the command-line option."""

    grid = create_dummy_grid()
    config = grid.ctx().config()

    L = config.get_double("water_latent_heat_fusion")
    rho = config.get_double("ice_density")

    # prescribe a heat flux that corresponds to a mass flux which is
    # an integer multiple of m / year so that we can easily specify it
    # using a command-line option
    M = PISM.convert(grid.ctx().unit_system(), 1, "m / year", "m / second")
    Q_default = config.get_double("ocean_sub_shelf_heat_flux_into_ice")
    Q = M * L * rho
    config.set_double("ocean_sub_shelf_heat_flux_into_ice", Q)

    # without the command-line option
    ocean_constant = PISM.OceanConstant(grid)
    ocean_constant.init()
    mass_flux_1 = PISM.model.createIceThicknessVec(grid)
    ocean_constant.shelf_base_mass_flux(mass_flux_1)

    # reset Q
    config.set_double("ocean_sub_shelf_heat_flux_into_ice", Q_default)

    # with the command-line option
    o = PISM.PETSc.Options()
    o.setValue("-shelf_base_melt_rate", 1.0)

    ocean_constant = PISM.OceanConstant(grid)
    ocean_constant.init()
    mass_flux_2 = PISM.model.createIceThicknessVec(grid)
    ocean_constant.shelf_base_mass_flux(mass_flux_2)

    import numpy as np
    with PISM.vec.Access(nocomm=[mass_flux_1, mass_flux_2]):
        assert np.fabs(mass_flux_1[0, 0] - M * rho) < 1e-16
        assert np.fabs(mass_flux_2[0, 0] - M * rho) < 1e-16

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

    def compare(pio):
        pio.open(basename + ".nc", PISM.PISM_READONLY)
        read_string = pio.get_att_text("PISM_GLOBAL", "string_attribute")
        read_text = pio.get_att_text("PISM_GLOBAL", "text_attribute")
        pio.close()

        # check that written and read strings are the same
        print "written string: '%s'" % attribute
        print "read string:    '%s'" % read_string
        print "read text:      '%s'" % read_text
        assert read_string == attribute
        assert read_text == attribute

    def netcdf3():
        # try reading this attribute
        pio = PISM.PIO(PISM.PETSc.COMM_WORLD, "netcdf3")

        print "\nTesting pism::io::NC3File::get_att_text_impl()..."
        compare(pio)

    def netcdf4_parallel():
        # try reading this attribute
        try:
            # try creating a netcdf4_parallel backend, stop the test
            # (without failing) if this fails -- PISM may not have
            # been compiled with parallel NetCDF-4.
            pio = PISM.PIO(PISM.PETSc.COMM_WORLD, "netcdf4_parallel")
        except:
            return

        print "\nTesting pism::io::NC4File::get_att_text_impl()..."
        compare(pio)

    setup()

    netcdf3()
    netcdf4_parallel()

    teardown()

netcdf_string_attribute_test()
