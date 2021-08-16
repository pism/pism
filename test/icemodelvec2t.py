# run this using "python -m nose icemodelvec2t.py -m test_name" to run only one of the
# tests
import unittest
import numpy
import os
import PISM
from PISM.testing import filename

ctx = PISM.Context()

# use the 360-day calendar so that the period of 1 year has a clear meaning
ctx.config.set_string("time.calendar", "360_day")
ctx.ctx.time().init_calendar("360_day")

# suppress all output
ctx.log.set_threshold(1)

def compare(v, value):
    with PISM.vec.Access(nocomm=v):
        numpy.testing.assert_almost_equal(v[0, 0], value)

class ForcingInput(unittest.TestCase):

    def setUp(self):
        "Prepare an input file with an interesting time-dependent field."
        self.filename = filename("input_")
        self.empty = filename("empty_")
        self.one_record = filename("one_record_")
        self.no_time = filename("no_time_")
        self.no_bounds = filename("no_time_bounds_")
        self.time_order = filename("time_order_")

        M = 3
        self.grid = PISM.IceGrid.Shallow(ctx.ctx, 1, 1, 0, 0, M, M, PISM.CELL_CORNER,
                                         PISM.NOT_PERIODIC)

        v = PISM.IceModelVec2S(self.grid, "v", PISM.WITHOUT_GHOSTS)

        units = "days since 1-1-1"
        N = 12
        dt = 30.0                 # days (floating point)
        t = numpy.r_[0:N] * dt
        self.tb = numpy.r_[0:N + 1] * dt
        self.f = numpy.sin(2.0 * numpy.pi * t / (N * dt))

        def write_data(filename, use_bounds=True, forward=True):
            bounds = PISM.VariableMetadata("time_bounds", ctx.unit_system)
            bounds.set_string("units", units)

            output = PISM.util.prepare_output(filename, append_time=False)
            output.write_attribute("time", "units", units)

            if use_bounds:
                PISM.define_time_bounds(bounds, "time", "nv", output, PISM.PISM_DOUBLE)
                output.write_attribute("time", "bounds", "time_bounds")

            if forward:
                order = range(N)
            else:
                order = range(N - 1, -1, -1)

            for k in order:
                PISM.append_time(output, "time", self.tb[k+1])

                if use_bounds:
                    PISM.write_time_bounds(output, bounds, k, (self.tb[k], self.tb[k + 1]))

                v.set(self.f[k])
                v.write(output)

        # regular forcing file
        write_data(self.filename, use_bounds=True)

        # file without time bounds
        write_data(self.no_bounds, use_bounds=False)

        # file with invalid time order
        write_data(self.time_order, forward=False)

        # empty file
        PISM.util.prepare_output(self.empty)

        # file with one record
        output = PISM.util.prepare_output(self.one_record)
        v.set(self.f[-1])
        v.write(output)

        # file without a time dimension
        v.set(self.f[-1])
        v.set_time_independent(True)
        v.dump(self.no_time)

    def tearDown(self):
        "Remove files created by setUp()"
        files = [self.filename, self.empty, self.one_record,
                 self.no_time, self.no_bounds, self.time_order]
        for f in files:
            os.remove(f)

    def forcing(self, filename, buffer_size=12, periodic=False):
        "Allocate and initialize forcing"
        input_file = PISM.File(ctx.com, self.filename, PISM.PISM_NETCDF3, PISM.PISM_READONLY)
        forcing = PISM.IceModelVec2T.ForcingField(self.grid, input_file, "v", "",
                                                  buffer_size, 52, periodic)
        forcing.metadata().set_string("long_name", "test field")

        forcing.init(filename, periodic)

        return forcing

    def test_missing_file(self):
        "Missing input file"
        try:
            forcing = self.forcing(self.empty)
            assert 0, "initialized with an empty file"
        except RuntimeError:
            pass

    def test_average(self):
        "Time average using average(t, dt)"
        forcing = self.forcing(self.filename)

        # second month
        N = 1
        t = self.tb[N] * 86400 + 1
        dt = forcing.max_timestep(t).value()
        forcing.update(t, dt)
        forcing.average(t, dt)

        compare(forcing, self.f[N])

        # average with a zero length
        forcing.average(t, 0.0)

        compare(forcing, self.f[N])

        # average() with only one record
        forcing = self.forcing(self.filename, buffer_size=1)
        forcing.update(0, 1)
        forcing.average(0, 1)

        compare(forcing, self.f[0])

    def test_extrapolation_left(self):
        "Extrapolation on the left"
        forcing = self.forcing(self.filename)

        t = self.tb[0] * 86400 - 1
        dt = self.tb[1] * 86400 + 1
        forcing.update(0, dt)
        forcing.interp(t)

        compare(forcing, self.f[0])

    def test_extrapolation_right(self):
        "Extrapolation on the right"
        forcing = self.forcing(self.filename)

        t = self.tb[-1] * 86400
        forcing.update(0, t)
        forcing.interp(t + 1)

        compare(forcing, self.f[-1])

    def test_interp_1(self):
        "Interpolation using interp(t)"
        forcing = self.forcing(self.filename)

        N = 1
        t = self.tb[N] * 86400 + 1
        forcing.update(0, t)
        forcing.interp(t)

        compare(forcing, self.f[N])

    def test_interp_2(self):
        "Interpolation using init_interpolation(ts) and interp(i, j)"
        forcing = self.forcing(self.filename)
        N = 12
        dt = 30.0                 # days (floating point)
        ts = numpy.r_[0:N] * dt + 0.5 * dt
        ts = [PISM.util.convert(T, "days", "seconds") for T in ts]

        forcing.update(0, self.tb[-1])
        forcing.init_interpolation(ts)

        with PISM.vec.Access(nocomm=forcing):
            numpy.testing.assert_almost_equal(forcing.interp(0, 0), self.f)

    def test_one_record(self):
        "Input file with only one time record"
        forcing = self.forcing(self.one_record)

        self.check_forcing(forcing, self.f[-1], 0, 1)

    def test_no_time_dimension(self):
        "Forcing without a time dimension"
        forcing = self.forcing(self.no_time)

        self.check_forcing(forcing, self.f[-1], 0, 1)

    def test_constant_field(self):
        "Field allocatted using IceModelVec2T::Constant()"
        f = 100.0
        forcing = PISM.IceModelVec2T.Constant(self.grid, "v", f)

        self.check_forcing(forcing, f, 0, 1)

    def check_forcing(self, forcing, value, t, dt):
        forcing.update(t, dt)
        forcing.average(t, dt)
        compare(forcing, value)

    def test_large_update_interval(self):
        "Test update() calls with an update interval past the last available time"
        forcing = self.forcing(self.filename)

        forcing.update(1e12, 1e12+1)
        forcing.init_interpolation([1e12])

        with PISM.vec.Access(nocomm=forcing):
            numpy.testing.assert_almost_equal(forcing.interp(0, 0), self.f[-1])

    def test_buffer_too_small(self):
        "Reading periodic data that does not fit in the buffer"
        # Note: IceModelVec2T::init() will throw RuntimeError if the buffer is too small
        # to hold all records of a periodic forcing field. This will never happen if this
        # IceModelVec2T was allocated using IceModelVec2T::ForcingField because it ensures
        # that the buffer is big enough.
        forcing = self.forcing(self.filename, buffer_size=2, periodic=False)
        try:
            forcing.init(self.filename, periodic=True)
            assert False, "Failed to stop because the buffer size is too small"
        except RuntimeError:
            pass

    def test_time_step_too_long(self):
        "update() call with a time step that is too long"

        N = 2
        forcing = self.forcing(self.filename, buffer_size=N)

        try:
            dt = (N + 1) * 30 * 86400
            forcing.update(0, dt)
            assert False, "Failed to catch a time step that is too long"
        except RuntimeError:
            pass

    def test_no_time_bounds(self):
        "Forcing without time bounds"
        try:
            forcing = self.forcing(self.no_bounds)
            assert False, "Loaded forcing without time bounds"
        except RuntimeError:
            pass

    def test_decreasing_time(self):
        "Invalid (decreasing) time"
        try:
            forcing = self.forcing(self.time_order)
            assert False, "Loaded forcing with a decreasing time dimension"
        except RuntimeError:
            pass

    def test_multiple_steps(self):
        "miltiple update() calls with a small buffer, discarding records"
        # In this test it is crucial that the time step (dt below) is long enough to use
        # more than one time record. This triggers the code that discards records that are
        # no longer needed but keeps the ones that are still of use, reducing the amount
        # of data we have to read from a file.

        forcing = self.forcing(self.filename, buffer_size=3)

        def check(month):
            t = self.tb[month] * 86400 + 1
            dt = 2 * 30 * 86400     # long-ish time step (2 months)
            forcing.update(t, dt)
            forcing.interp(t)

            compare(forcing, self.f[month])

        # second month
        check(1)
        # fourth month
        check(3)

    def test_max_timestep(self):
        "Maximum time step"
        forcing = self.forcing(self.filename, buffer_size=1)

        # time bounds in seconds
        tb = self.tb * 86400

        assert forcing.max_timestep(0).value() == tb[1] - tb[0]
        assert forcing.max_timestep(tb[-1]).infinite()
        assert forcing.max_timestep(tb[1] - 0.5).value() == tb[2] - tb[1]
        assert forcing.max_timestep(tb[-1] - 0.5).infinite()

    def test_periodic(self):
        "Using periodic data"
        forcing = self.forcing(self.filename, periodic=True)

        # a year and a half
        t = numpy.arange(18 + 1) * 30 * 86400.0

        forcing.update(0, t[-1])

        forcing.init_interpolation(t)

        with PISM.vec.Access(nocomm=forcing):
            a = numpy.array(forcing.interp(0, 0))
        b = numpy.r_[self.f, self.f[0:(6 + 1)]]

        numpy.testing.assert_almost_equal(a, b)
