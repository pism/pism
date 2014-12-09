"""This file contains very superficial tests of the PISM Python
wrappers. The goal is to be able to detect changes in the API
(function signatures, etc), not to test correctness.

Use with nose (https://pypi.python.org/pypi/nose/) and coverage.py
(https://pypi.python.org/pypi/coverage)

Run this to get a coverage report:

nosetests --with-coverage --cover-branches --cover-html --cover-package=PISM test/nosetests.py
"""

import PISM

def create_dummy_grid():
    "Create a dummy grid"
    ctx = PISM.Context()
    grid = ctx.newgrid()
    PISM.model.initShallowGrid(grid, 1e5, 1e5, 100, 100, PISM.NOT_PERIODIC)
    return grid

def context_test():
    "Test creating a new PISM context"
    ctx = PISM.Context()
    config = ctx.config

def context_missing_attribute_test():
    "Test the handling of missing attributes"
    ctx = PISM.Context()
    try:
        config = ctx.foo        # there is no "foo", this should fail
        return False
    except AttributeError:
        return True

def init_config_twice_test():
    "Test initializing the config twice"
    ctx = PISM.Context()
    config = ctx.config         # this should initialize config

    try:
        ctx.init_config()       # initialize again; this should fail
        return False
    except RuntimeError:
        return True

def create_grid_test():
    "Test the creation of the IceGrid object"
    create_dummy_grid()

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

def vec_access_test():
    "Test the PISM.vec.Access class and IceGrid::points, points_with_ghosts, coords"
    grid = create_dummy_grid()

    vec_scalar = PISM.vec.randVectorS(grid, 1.0)
    vec_scalar_ghosted = PISM.vec.randVectorS(grid, 1.0, 2)

    with PISM.vec.Access(comm=[vec_scalar_ghosted], nocomm=vec_scalar):
        for (i,j) in grid.points_with_ghosts():
            pass

    with PISM.vec.Access(comm=vec_scalar_ghosted, nocomm=[vec_scalar]):
        for (i,j) in grid.points():
            # do something
            pass

        for (i,j,x,y) in grid.coords():
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

    tz2 = PISM.vec.ToProcZero(grid,dof=2,dim=2)
    array_vector_0 = tz2.communicate(vec_vector)

    try:
        tz3 = PISM.vec.ToProcZero(grid,dof=2,dim=3)
        return False
    except NotImplementedError:
        # 3D fields are not supported (yet)
        pass

def create_modeldata_test():
    "Test creating the ModelData class"
    grid = create_dummy_grid()
    md = PISM.model.ModelData(grid)

def grid_from_file_test():
    "Intiialize a grid from a file"
    grid = create_dummy_grid()

    enthalpy = PISM.model.createEnthalpyVec(grid)
    enthalpy.set(80e3)

    output_file = "test_grid_from_file.nc"
    pio = PISM.PIO(grid,"netcdf3")
    pio.open(output_file, PISM.PISM_READWRITE_MOVE)
    pio.def_time(grid.config.get_string("time_dimension_name"),
                 grid.config.get_string("calendar"), grid.time.units_string())
    pio.append_time(grid.config.get_string("time_dimension_name"),grid.time.current())
    pio.close()

    enthalpy.write(output_file)

    grid2 = PISM.Context().newgrid()

    PISM.model.initGridFromFile(grid, output_file, PISM.NOT_PERIODIC)

def create_special_vecs_test():
    "Test helpers used to create standard PISM fields"
    grid = create_dummy_grid()

    usurf = PISM.model.createIceSurfaceVec(grid)

    thk = PISM.model.createIceThicknessVec(grid)

    usurfstore = PISM.model.createIceSurfaceStoreVec(grid)

    thkstore = PISM.model.createIceThicknessStoreVec(grid)

    bed = PISM.model.createBedrockElevationVec(grid)

    tauc = PISM.model.createYieldStressVec(grid)

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

    # test var_cmp
    PISM.model.var_cmp(lon, lat)
    PISM.model.var_cmp(lat, lon)
    PISM.model.var_cmp(lon, lon)

    # test ModelVecs.add()
    modeldata = PISM.model.ModelData(grid)
    vecs = modeldata.vecs

    vecs.add(mask)

    return True

def options_test():
    "Test command-line option handling"
    ctx = PISM.Context()

    for o in PISM.OptionsGroup(ctx.com, "", "Testing options"):
        M = PISM.optionsInt("-M", "description", default=100)
        S = PISM.optionsString("-S", "description", default="string")
        R = PISM.optionsReal("-R", "description", default=1.5)
        B = PISM.optionsFlag("-B", "description", default=False)
        IA = PISM.optionsIntArray("-IA", "description", default=[1,2])
        RA = PISM.optionsRealArray("-RA", "description", default=[2,3])
        SA = PISM.optionsStringArray("-SA", "description", default="one,two")

    for o in PISM.OptionsGroup(None, "", "Trying None for the COMM"):
        M = PISM.optionsList(ctx.com, "-L", "description", opt_set="one,two", default="one")
        # without the default
        SA = PISM.optionsStringArray("-SA", "description")

def modelvecs_test():
    "Test the ModelVecs class"

    ctx = PISM.Context()
    grid = ctx.newgrid()
    PISM.model.initGrid(grid, 1e5, 1e5, 1000.0, 100, 100, 11, PISM.NOT_PERIODIC)

    mask = PISM.model.createIceMaskVec(grid)
    mask.set(PISM.MASK_GROUNDED)

    modeldata = PISM.model.ModelData(grid)
    vecs = modeldata.vecs

    vecs.add(mask, "ice_mask", writing=True)

    try:
        vecs.add(mask, "ice_mask")
        return False
    except RuntimeError:
        # should fail: mask was added already
        pass

    # get a field:
    vecs.get("ice_mask")

    try:
        vecs.get("invalid")
        return False
    except AttributeError:
        # should fail
        pass

    # test __repr__
    print vecs

    # test rename()
    vecs.rename("ice_mask", "mask")

    # remove() a vec not marked for writing
    vecs.remove("mask")

    vecs.add(mask, "mask")

    # test setPISMVarsName
    vecs.setPISMVarsName("mask", "pism_mask")

    try:
        vecs.setPISMVarsName("invalid", "new_name")
        return False
    except RuntimeError:
        # should fail: cannot set the name of a variable that is not
        # there
        pass

    # test asPISMVars()
    vecs.asPISMVars()

    # test has()
    vecs.has("thickness")

    # test markForWriting
    vecs.markForWriting("mask")
    vecs.markForWriting(mask)

    # test write()
    output_file = "test_ModelVecs.nc"
    pio = PISM.PIO(grid,"netcdf3")
    pio.open(output_file, PISM.PISM_READWRITE_MOVE)
    pio.def_time(grid.config.get_string("time_dimension_name"),
                 grid.config.get_string("calendar"), grid.time.units_string())
    pio.append_time(grid.config.get_string("time_dimension_name"),grid.time.current())
    pio.close()

    vecs.write(output_file)

    # test writeall()
    vecs.writeall(output_file)

    # test remove()
    vecs.remove("mask")


def sia_test():
    "Test the PISM.sia module"
    ctx = PISM.Context()
    grid = ctx.newgrid()
    PISM.model.initGrid(grid, 1e5, 1e5, 1000.0, 100, 100, 11, PISM.NOT_PERIODIC)

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
    enthalpy.set(enthalpyconverter.getEnth(270.0, 0.0, 0.0))

    modeldata = PISM.model.ModelData(grid)
    modeldata.setPhysics(enthalpyconverter)

    vecs = modeldata.vecs;

    fields = {"thickness": thk,
              "surface": surface,
              "ice_mask": mask,
              "bed": bed,
              "enthalpy": enthalpy}

    for key in fields.keys():
        vecs.add(fields[key], key)

    vel_sia = PISM.sia.computeSIASurfaceVelocities(modeldata)

def util_test():
    "Test the PISM.util module"
    grid = create_dummy_grid()

    output_file = "test_pism_util.nc"
    pio = PISM.PIO(grid,"netcdf3")
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
    pio = PISM.PIO(grid, "netcdf3")

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
    c.write("other_log.nc", "other_log") # non-default arguments

def column_interpolation_test(plot=False):
    """Test ColumnInterpolation by interpolating from the coarse grid to the
    fine grid and back."""
    import numpy as np
    import pylab as plt

    Lz = 1000.0
    Mz = 41

    def z_quadratic(Mz, Lz):
        result = np.zeros(Mz)
        z_lambda = 4.0
        for k in xrange(Mz-1):
            zeta = float(k) / (Mz - 1)
            result[k] = Lz * ((zeta / z_lambda) * (1.0 + (z_lambda - 1.0) * zeta));
        result[Mz-1] = Lz
        return result

    def test_quadratic_interp():
        z_coarse = z_quadratic(Mz, Lz)
        f_coarse = (z_coarse / Lz)**2

        print "Testing quadratic interpolation"
        return test_interp(z_coarse, f_coarse)

    def test_linear_interp():
        z_coarse = np.linspace(0, Lz, Mz)
        f_coarse = (z_coarse / Lz)**2

        print "Testing linear interpolation"
        return test_interp(z_coarse, f_coarse)

    def test_interp(z, f):
        interp = PISM.ColumnInterpolation(z)

        z_fine = interp.z_fine()

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
            plt.grid(True)

        if plot:
            plot()

        delta = np.linalg.norm(f - f_roundtrip, ord=1)
        delta_numpy = np.linalg.norm(f_fine - f_fine_numpy, ord=1)
        print "norm1(fine_to_coarse(coarse_to_fine(f)) - f) = %f" % delta
        print "norm1(PISM - NumPy) = %f" % delta_numpy

        return delta,delta_numpy

    linear_delta,linear_delta_numpy = test_linear_interp()

    quadratic_delta,_ = test_quadratic_interp()

    if plot:
        plt.show()

    if (linear_delta > 1e-12 or
        linear_delta_numpy > 1e-12 or
        quadratic_delta > 1e-3):
        return False
    return True
