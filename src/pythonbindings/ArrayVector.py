def __getitem__(self, *args):
    i, j = args[0]
    return self.getitem(i, j)

def __setitem__(self, *args):
    if(len(args) == 2):
        i, j = args[0]
        value = args[1]
        if(isinstance(value, list) and len(value) == 2):
            u, v = value
            self.setitem(i, j, u, v)
        else:
            self.setitem(i, j, value)
    else:
        raise ValueError("__setitem__ requires 2 arguments; received %d" % len(args))

def to_xarray(self, **kwargs):
    "Return an xarray.Dataset (a copy) containing data from this vector field (on rank 0)."
    tmp = self.allocate_proc0_copy()
    self.put_on_proc0(tmp)
    import numpy
    import xarray
    import cftime

    if self.grid().ctx().rank() == 0:
        data = numpy.array(tmp.get()).reshape(1, *self.shape())
        grid = self.grid()
        time = grid.ctx().time()
        calendar = time.calendar()
        units = time.units_string()
        date = cftime.num2date(time.current(), units, calendar)
        spatial_coords = self.spatial_coords
        das = []
        for k in range(data.shape[-1]):
            metadata = self.metadata(k)
            name = metadata.get_string("short_name")
            if metadata.n_spatial_dimensions() > 2:
                dims = ["time", "y", "x", "z"]
            else:
                dims = ["time", "y", "x"]
            
            attrs = self.attrs
            # Spatial coordinates
            coords = {k: ([k], numpy.array(grid[k]), spatial_coords[k]) for k in spatial_coords.keys()}
            # Temporal coordinate
            coords["time"] =  (["time"], xarray.cftime_range(date, periods=1), {})
            das.append(xarray.DataArray(data[...,k], coords=coords, dims=dims, attrs=attrs, name=name, **kwargs))
        return xarray.merge(das)
    else:
        return None
