import PISM
import PISM.testing
import os
from unittest import TestCase, SkipTest

# always available
backends = [PISM.PISM_NETCDF3]

if PISM.Pism_USE_PARALLEL_NETCDF4:
    backends += [PISM.PISM_NETCDF4_PARALLEL]

if PISM.Pism_USE_PNETCDF:
    backends += [PISM.PISM_PNETCDF]

if PISM.Pism_USE_PIO:
    # assume that ParallelIO was built with all possible libraries
    backends += [PISM.PISM_PIO_PNETCDF, PISM.PISM_PIO_NETCDF,
                 PISM.PISM_PIO_NETCDF4C, PISM.PISM_PIO_NETCDF4P]

ctx = PISM.Context().ctx

backend_names = {PISM.PISM_NETCDF3 : "netcdf3",
                 PISM.PISM_NETCDF4_PARALLEL : "netcdf4_parallel",
                 PISM.PISM_PNETCDF : "pnetcdf",
                 PISM.PISM_PIO_NETCDF : "pio_netcdf",
                 PISM.PISM_PIO_NETCDF4P : "pio_netcdf4p",
                 PISM.PISM_PIO_NETCDF4C : "pio_netcdf4c",
                 PISM.PISM_PIO_PNETCDF: "pio_pnetcdf"}

def fail(backend):
    assert False, "test failed (backend = {})".format(backend_names[backend])

def test_string_to_backend():
    "PISM.string_to_backend()"

    for backend, name in backend_names.items():
        assert PISM.string_to_backend(name) == backend

    try:
        PISM.string_to_backend("invalid")
        return False
    except RuntimeError:
        pass

class File(TestCase):

    def test_empty_filename(self):
        for backend in backends:
            try:
                f = PISM.File(ctx.com(), "", backend, PISM.PISM_READONLY,
                              ctx.pio_iosys_id())
                fail(backend)
            except RuntimeError:
                pass

    def test_missing_file(self):
        for backend in backends:
            try:
                f = PISM.File(ctx.com(), "missing_file.nc", backend, PISM.PISM_READONLY,
                              ctx.pio_iosys_id())
                fail(backend)
            except RuntimeError:
                pass

    def test_backend_guessing(self):
        "File(..., PISM_GUESS, ...)"

        f = PISM.File(ctx.com(), self.file_with_time, PISM.PISM_GUESS, PISM.PISM_READONLY,
                      ctx.pio_iosys_id())
        assert f.nrecords() == 1

    def test_create_clobber(self):
        "File(..., PISM_READWRITE_CLOBBER)"
        try:
            for backend in backends:
                f = PISM.File(ctx.com(), "test_filename.nc", backend, PISM.PISM_READWRITE_CLOBBER,
                              ctx.pio_iosys_id())
        finally:
            os.remove("test_filename.nc")

    def test_create_move(self):
        "File(..., PISM_READWRITE_MOVE)"
        try:
            for backend in backends:
                f = PISM.File(ctx.com(), "test_filename.nc", backend, PISM.PISM_READWRITE_MOVE,
                              ctx.pio_iosys_id())
        finally:
            os.remove("test_filename.nc")
            try:
                os.remove("test_filename.nc~")
            except:
                pass

    def test_com(self):
        "File.com()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())
            assert f.com() == ctx.com()
            f.close()

    def test_close(self):
        "File.close()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())
            f.close()

            try:
                f.close()
                # closing twice is an error
                fail(backend)
            except RuntimeError:
                pass

    def test_redef(self):
        "File.redef()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READWRITE,
                          ctx.pio_iosys_id())
            f.redef()
            f.close()

    def test_enddef(self):
        "File.enddef()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READWRITE,
                          ctx.pio_iosys_id())
            f.enddef()
            f.close()

    def test_sync(self):
        "File.sync()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READWRITE,
                          ctx.pio_iosys_id())
            f.sync()
            f.close()

    def test_name(self):
        "File.name()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())
            assert f.name() == self.file_with_time
            f.close()

    def test_nrecords(self):
        "File.nrecords()"
        for F in self.files:
            for backend in backends:
                f = PISM.File(ctx.com(), F, backend, PISM.PISM_READONLY,
                              ctx.pio_iosys_id())
                assert f.nrecords() == 1
                f.close()

    def test_nrecords_variable(self):
        "File.nrecords(variable)"
        for F in [self.file_with_time, self.file_without_time]:
            for backend in backends:
                f = PISM.File(ctx.com(), F, backend, PISM.PISM_READONLY,
                              ctx.pio_iosys_id())
                assert f.nrecords("v", "standard_name", ctx.unit_system()) == 1
                # found using the standard name
                assert f.nrecords("w", "standard_name", ctx.unit_system()) == 1
                assert f.nrecords("missing", "", ctx.unit_system()) == 0
                f.close()

    def test_nvariables(self):
        "File.nvariables()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())
            assert f.nvariables() == 4 # time, x, y, v
            f.close()

    def test_nattributes(self):
        "File.nattributes()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())
            assert f.nattributes("time") == 4 # units, axis, calendar, long_name
            assert f.nattributes("x") == 5 # units, axis, long_name, standard_name, spacing_meters
            assert f.nattributes("PISM_GLOBAL") == 2

            f.close()

    def test_define_dimension(self):
        "File.define_dimension()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READWRITE,
                          ctx.pio_iosys_id())
            f.define_dimension("dim_{}".format(backend), 10 + backend)
            f.close()

    def test_dimension_length(self):
        "File.dimension_length()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())
            assert f.dimension_length("time") == 1
            assert f.dimension_length("x") == 3
            assert f.dimension_length("y") == 5
            assert f.dimension_length("z") == 0
            f.close()

    def test_dimensions(self):
        "File.dimensions()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READWRITE,
                          ctx.pio_iosys_id())
            assert f.dimensions("v") == ("time", "y", "x")

            variable_name = "scalar_variable_{}".format(backend)
            f.define_variable(variable_name, PISM.PISM_BYTE, [])
            assert f.dimensions(variable_name) == ()
            f.close()

    def test_dimension_exists(self):
        "File.dimension_exists()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())
            assert f.dimension_exists("x")
            assert not f.dimension_exists("z")
            f.close()

    def test_dimension_type(self):
        "File.dimension_type()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())
            f.dimension_type("time", ctx.unit_system()) == PISM.T_AXIS
            f.dimension_type("x", ctx.unit_system()) == PISM.X_AXIS
            f.dimension_type("y", ctx.unit_system()) == PISM.Y_AXIS

            try:
                f.dimension_type("z", ctx.unit_system())
                fail(backend)
            except RuntimeError:
                pass

            f.close()

            f = PISM.File(ctx.com(), self.file_dim_types, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())

            def check(names, axis_type):
                for c in names:
                    assert f.dimension_type(c, ctx.unit_system()) == axis_type

            check(self.x_names + self.strange_x_names, PISM.X_AXIS)
            check(self.y_names + self.strange_y_names, PISM.Y_AXIS)
            check(self.z_names + self.strange_z_names, PISM.Z_AXIS)
            check(self.t_names + self.strange_t_names, PISM.T_AXIS)

            assert f.dimension_type("unknown_axis", ctx.unit_system()) == PISM.UNKNOWN_AXIS

            f.close()

    def test_read_dimension(self):
        "File.read_dimension()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())
            assert f.read_dimension("x") == (-10000.0, 0.0, 10000.0)

            try:
                f.read_dimension("z")
                fail(backend)
            except RuntimeError:
                pass

            f.close()

    def test_variable_name(self):
        "File.variable_name()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())
            assert f.variable_name(0) == "time"

            # invalid variable index
            try:
                f.variable_name(1000)
                fail(backend)
            except RuntimeError:
                pass

            f.close()

    def test_define_variable(self):
        "File.define_variable()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READWRITE,
                          ctx.pio_iosys_id())
            f.define_variable("var_{}".format(backend), PISM.PISM_DOUBLE, ["y", "x"])
            # defining an existing variable should fail
            try:
                f.define_variable("var_{}".format(backend), PISM.PISM_DOUBLE, ["y", "x"])
                fail(backend)
            except RuntimeError:
                pass

            # defining a variable depending on a non-existent dimension should fail gracefully
            try:
                f.define_variable("var_{}_1".format(backend), PISM.PISM_DOUBLE, ["y", "x", "z"])
                fail(backend)
            except RuntimeError:
                pass
            f.close()

    def test_find_variable(self):
        "File.find_variable(short_name)"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())
            assert f.find_variable("time")
            assert not f.find_variable("z")
            f.close()

    def test_find_variable(self):
        "File.find_variable(short_name, standard_name)"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())
            assert f.find_variable("v", "standard_name").exists
            assert f.find_variable("v", "standard_name").found_using_standard_name
            assert f.find_variable("other_name", "standard_name").found_using_standard_name
            assert f.find_variable("other_name", "standard_name").name == "v"
            assert not f.find_variable("v", "").found_using_standard_name
            assert f.find_variable("missing", "other_standard_name").exists == False
            assert f.find_variable("missing", "other_standard_name").name == ""
            f.close()

            f = PISM.File(ctx.com(), self.file_inconsistent, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())

            try:
                f.find_variable("v", "standard_name")
                fail(backend)
            except RuntimeError:
                pass

    def test_read_variable(self):
        "File.read_variable()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())
            f.read_variable("x", [1], [1]) == [0.0]

            # start and count of different lengths
            if PISM.Pism_DEBUG:
                try:
                    f.read_variable("v", [1, 1], [1, 1, 1])
                    fail(backend)
                except RuntimeError:
                    pass

            f.close()

    def test_read_variable_transposed(self):
        raise SkipTest("disable this in Python bindings")

    def test_write_variable(self):
        "File.write_variable()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_without_time, backend, PISM.PISM_READWRITE,
                          ctx.pio_iosys_id())
            f.write_variable("v", [1, 1], [1, 1], [100.0])
            f.sync()
            assert f.read_variable("v", [1, 1], [1, 1]) == (100.0,)

            # start and count of different lengths
            if PISM.Pism_DEBUG:
                try:
                    f.write_variable("v", [1, 1], [1, 1, 1], [200.0])
                    fail(backend)
                except RuntimeError:
                    pass

            f.close()

    def test_write_distributed_array(self):
        raise SkipTest("disable this in Python bindings")

    def test_remove_attribute(self):
        "File.remove_attribute()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READWRITE,
                          ctx.pio_iosys_id())
            assert f.attribute_type("time", "units") == PISM.PISM_CHAR
            f.remove_attribute("time", "units")
            assert f.attribute_type("time", "units") == PISM.PISM_NAT
            f.write_attribute("time", "units", "seconds since 1-1-1")
            assert f.attribute_type("time", "units") == PISM.PISM_CHAR
            f.close()

    def test_attribute_name(self):
        "File.attribute_name()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())
            assert f.attribute_name("time", 0) == "units"
            assert f.attribute_name("PISM_GLOBAL", 0) == "global_text_attr"
            f.close()

    def test_attribute_type(self):
        "File.attribute_type()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())
            f.attribute_type("x", "units") == PISM.PISM_CHAR
            f.attribute_type("x", "spacing_meters") == PISM.PISM_DOUBLE
            f.attribute_type("x", "missing") == PISM.PISM_NAT
            f.attribute_type("PISM_GLOBAL", "global_text_att") == PISM.PISM_CHAR
            f.close()

    def test_write_attribute_number(self):
        "File.write_attribute(number)"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READWRITE,
                          ctx.pio_iosys_id())
            f.write_attribute("v", "new_attribute", PISM.PISM_DOUBLE, (1.0, 2.0))
            assert f.read_double_attribute("v", "new_attribute") == (1.0, 2.0)

            f.write_attribute("PISM_GLOBAL", "new_attribute", PISM.PISM_DOUBLE, (1.0, 2.0))
            assert f.read_double_attribute("PISM_GLOBAL", "new_attribute") == (1.0, 2.0)

            f.close()

    def test_write_attribute_string(self):
        "File.write_attribute(string)"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READWRITE,
                          ctx.pio_iosys_id())
            f.write_attribute("v", "new_attribute", "test string")
            assert f.read_text_attribute("v", "new_attribute") == "test string"

            f.write_attribute("PISM_GLOBAL", "new_global_attr", "test_global")
            assert f.read_text_attribute("PISM_GLOBAL", "new_global_attr") == "test_global"
            f.close()

    def test_read_double_attribute(self):
        "File.read_double_attribute()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())
            assert f.read_double_attribute("x", "spacing_meters") == (10000.0,)
            assert f.read_double_attribute("x", "missing") == ()
            # type mismatch: fail with a helpful message
            try:
                f.read_double_attribute("x", "units")
                fail(backend)
            except RuntimeError:
                pass

            assert f.read_double_attribute("PISM_GLOBAL", "global_double_attr") == (12.0,)

            f.close()

    def test_read_text_attribute(self):
        "File.read_text_attribute()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READONLY,
                          ctx.pio_iosys_id())
            assert f.read_text_attribute("x", "units") == "m"
            assert f.read_text_attribute("x", "missing") == ""

            assert f.read_text_attribute("PISM_GLOBAL", "global_text_attr") == "test_global"

            # throw an error in case of a type mismatch
            try:
                f.read_text_attribute("x", "spacing_meters")
                fail(backend)
            except RuntimeError:
                pass

            f.close()

    def test_append_history(self):
        "File.read_text_attribute()"
        for backend in backends:
            f = PISM.File(ctx.com(), self.file_with_time, backend, PISM.PISM_READWRITE,
                          ctx.pio_iosys_id())
            try:
                f.remove_attribute("PISM_GLOBAL", "history")
            except:
                pass
            f.append_history("one")
            f.append_history("two")
            assert f.read_text_attribute("PISM_GLOBAL", "history") == "twoone"
            f.close()

    def setUp(self):
        self.file_with_time = "test_file_with_time.nc"
        self.file_without_time = "test_file_without_time.nc"
        self.file_inconsistent = "test_file_inconsistent.nc"
        self.file_dim_types = "test_file_dim_types.nc"

        self.files = [self.file_with_time, self.file_without_time, self.file_inconsistent, self.file_dim_types]

        grid = PISM.testing.shallow_grid()
        vec = PISM.Scalar(grid, "v")
        vec.metadata(0).long_name("dummy variable for testing").units("kelvin").output_units("degree_Celsius").standard_name("standard_name")
        vec.set(1.0)
        vec.metadata().set_time_independent(True)
        vec.dump(self.file_without_time)

        # file with two variables with the same standard name
        vec.dump(self.file_inconsistent)
        vec.metadata(0).set_name("w")
        vec.write(self.file_inconsistent)
        vec.metadata(0).set_name("v")

        vec.set(2.0)
        vec.metadata().set_time_independent(False)
        vec.dump(self.file_with_time)

        f = PISM.File(ctx.com(), self.file_with_time, PISM.PISM_NETCDF3, PISM.PISM_READWRITE)
        f.write_attribute("PISM_GLOBAL", "global_text_attr", "test_global")
        f.write_attribute("PISM_GLOBAL", "global_double_attr", PISM.PISM_DOUBLE, [12.0])
        f.close()

        f = PISM.File(ctx.com(), self.file_dim_types, PISM.PISM_NETCDF3, PISM.PISM_READWRITE_CLOBBER)
        self.x_names = ["x", "X", "x1", "X1"]
        strange_x_attrs = {"coord_x_1" : ("axis", "x"),
                           "coord_x_2" : ("axis", "X")}
        self.strange_x_names = list(strange_x_attrs.keys())

        self.y_names = ["y", "Y", "y1", "Y1"]
        strange_y_attrs = {"coord_y_1" : ("axis", "y"),
                           "coord_y_2" : ("axis", "Y")}
        self.strange_y_names = list(strange_y_attrs.keys())

        self.z_names = ["z", "Z", "z1", "Z1"]
        strange_z_attrs = {"coord_z_1" : ("axis", "z"),
                           "coord_z_2" : ("axis", "Z")}
        self.strange_z_names = list(strange_z_attrs.keys())

        self.t_names = ["t", "T", "time", "t0", "T0"]
        strange_t_attrs = {"coord_t_1" : ("units", "years"),
                           "coord_t_2" : ("standard_name", "time"),
                           "coord_t_3" : ("axis", "T"),
                           "coord_t_4" : ("axis", "t")}
        self.strange_t_names = list(strange_t_attrs.keys())

        def create(names, strange_names, attrs):
            for v in names + strange_names:
                f.define_variable(v, PISM.PISM_DOUBLE, [])

            for c, (a, v) in attrs.items():
                f.write_attribute(c, a, v)

        create(self.x_names, self.strange_x_names, strange_x_attrs)
        create(self.y_names, self.strange_y_names, strange_y_attrs)
        create(self.z_names, self.strange_z_names, strange_z_attrs)
        create(self.t_names, self.strange_t_names, strange_t_attrs)

        f.define_variable("unknown_axis", PISM.PISM_DOUBLE, [])

    def tearDown(self):
        for f in self.files:
            os.remove(f)
            pass

class StringAttribute(TestCase):
    "Test reading a NetCDF-4 string attribute."

    def test_read_string_attribute(self):
        "File.read_text_attribute() (string)"

        backends = [PISM.PISM_NETCDF3]

        if PISM.Pism_USE_PARALLEL_NETCDF4:
            backends += [PISM.PISM_NETCDF4_PARALLEL]

        for backend in backends:
            f = PISM.File(ctx.com(), self.basename + ".nc", backend, PISM.PISM_READONLY)
            assert self.attribute == f.read_text_attribute("PISM_GLOBAL", "string_attribute")
            assert self.attribute == f.read_text_attribute("PISM_GLOBAL", "text_attribute")
            # multi-valued string attributes are turned into comma-separated lists
            assert "{0},{0}".format(self.attribute) == f.read_text_attribute("PISM_GLOBAL",
                                                                             "string_multi_value")
            f.close()

    def setUp(self):
        self.basename = "string_attribute_test"
        self.attribute = "string attribute"

        cdl = """
netcdf string_attribute_test {{
  string :string_attribute = "{0}" ;
  string :string_multi_value = "{0}", "{0}" ;
  :text_attribute = "{0}" ;
}}
""".format(self.attribute)
        with open(self.basename + ".cdl", "w") as f:
            f.write(cdl)

        os.system("ncgen -4 %s.cdl" % self.basename)

    def tearDown(self):
        os.remove(self.basename + ".nc")
        os.remove(self.basename + ".cdl")
