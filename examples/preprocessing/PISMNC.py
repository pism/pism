#!/usr/bin/env python

try:
    import netCDF4 as netCDF
except:
    print("netCDF4 is not installed!")
    sys.exit(1)


class PISMDataset(netCDF.Dataset):

    def create_time(self, use_bounds=False, length=None, units=None):
        self.createDimension('time', size=length)
        t_var = self.createVariable('time', 'f8', ('time',))

        t_var.axis = "T"
        t_var.long_name = "time"
        if not units:
            t_var.units = "seconds since 1-1-1"  # just a default
        else:
            t_var.units = units

        if use_bounds:
            self.createDimension('n_bounds', 2)
            self.createVariable("time_bounds", 'f8', ('time', 'n_bounds'))
            t_var.bounds = "time_bounds"

    def create_dimensions(self, x, y, time_dependent=False, use_time_bounds=False):
        """
        Create PISM-compatible dimensions in a NetCDF file.
        """

        if time_dependent and not 'time' in list(self.variables.keys()):
            self.create_time(use_time_bounds)

        self.createDimension('x', x.size)
        self.createDimension('y', y.size)

        x_var = self.createVariable('x', 'f8', ('x',))
        x_var[:] = x

        y_var = self.createVariable('y', 'f8', ('y',))
        y_var[:] = y

        x_var.axis = "X"
        x_var.long_name = "X-coordinate in Cartesian system"
        x_var.units = "m"
        x_var.standard_name = "projection_x_coordinate"

        y_var.axis = "Y"
        y_var.long_name = "Y-coordinate in Cartesian system"
        y_var.units = "m"
        y_var.standard_name = "projection_y_coordinate"

        self.sync()

    def append_time(self, value, bounds=None):
        if 'time' in list(self.dimensions.keys()):
            time = self.variables['time']
            N = time.size
            time[N] = value

            if bounds:
                self.variables['time_bounds'][N, :] = bounds

    def set_attrs(self, var_name, attrs):
        """attrs should be a list of (name, value) tuples."""
        if not attrs:
            return

        for (name, value) in attrs.items():
            if name == "_FillValue":
                continue
            setattr(self.variables[var_name], name, value)

    def define_2d_field(self, var_name, time_dependent=False, dims=None, nc_type='f8', attrs=None):
        """
        time_dependent: boolean

        dims: an optional list of dimension names. use this to override the
              default order ('time', 'y', 'x')

        attrs: a dictionary of attributes
        """
        if not dims:
            if time_dependent:
                dims = ('time', 'y', 'x')
            else:
                dims = ('y', 'x')

        try:
            var = self.variables[var_name]
        except:
            if attrs is not None and '_FillValue' in list(attrs.keys()):
                var = self.createVariable(var_name, nc_type, dims,
                                          fill_value=attrs['_FillValue'])
            else:
                var = self.createVariable(var_name, nc_type, dims)

            self.set_attrs(var_name, attrs)

        return var

    def define_timeseries(self, var_name, attrs=None):
        try:
            if attrs is not None and '_FillValue' in list(attrs.keys()):
                var = self.createVariable(var_name, 'f8', ('time',),
                                          fill_value=attrs['_FillValue'])
            else:
                var = self.createVariable(var_name, 'f8', ('time',))
        except:
            var = self.variables[var_name]

        self.set_attrs(var_name, attrs)

        return var

    def write(self, var_name, data, time_dependent=False, attrs=None):
        """
        Write time-series or a 2D field to a file.
        """

        if data.ndim == 1:
            return self.write_timeseries(var_name, data, attrs=attrs)
        elif data.ndim == 2:
            return self.write_2d_field(var_name, data, time_dependent, attrs=attrs)
        else:
            return None

    def write_2d_field(self, var_name, data, time_dependent=False, attrs=None):
        """
        Write a 2D numpy array to a file in a format PISM can read.
        """

        var = self.define_2d_field(var_name, time_dependent, attrs=attrs)

        if time_dependent:
            last_record = self.variables['time'].size - 1
            var[last_record, :, :] = data
        else:
            var[:] = data

        return var

    def write_timeseries(self, var_name, data, attrs=None):
        """Write a 1D (time-series) array to a file."""

        var = self.define_timeseries(var_name, attrs=attrs)
        var[:] = data

        return var


if __name__ == "__main__":
    # produce a NetCDF file for testing
    from numpy import linspace, meshgrid

    nc = PISMDataset("foo.nc", 'w')

    x = linspace(-100, 100, 101)
    y = linspace(-100, 100, 201)

    xx, yy = meshgrid(x, y)

    nc.create_dimensions(x, y, time_dependent=True, use_time_bounds=True)

    nc.define_2d_field("xx", time_dependent=True,
                       attrs={"long_name": "xx",
                              "comment": "test variable",
                              "valid_range": (-200.0, 200.0)})

    for t in [0, 1, 2, 3]:
        nc.append_time(t, (t - 1, t))

        nc.write("xx", xx + t, time_dependent=True)
        nc.write("yy", yy + 2 * t, time_dependent=True)

    nc.close()
