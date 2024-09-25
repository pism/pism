# code extending the Array class


def regrid(self, filename, critical=False, default_value=0.0):
    if critical == True:
        self._regrid(filename, Default.Nil())
    else:
        self._regrid(filename, Default(default_value))

def to_numpy(self):
    "Return a NumPy array (a copy) containing data from this field (on rank 0)."
    tmp = self.allocate_proc0_copy()
    self.put_on_proc0(tmp)
    import numpy

    if self.grid().ctx().rank() == 0:
        return numpy.array(tmp.get()).reshape(self.shape())
    else:
        return None

def local_part(self):
    """NumPy array containing the local (sub-domain) part, ghosts and all.

    This is not a copy. Modifications have immediate effect.
    """

    padding = 2 * self.stencil_width()
    shape = self.shape()

    xm = self.grid().xm() + padding
    ym = self.grid().ym() + padding

    if len(shape) == 2:
        shape = (ym, xm)
    else:
        shape = (ym, xm, shape[2])

    return self.vec().get().array.reshape(shape)

@property
def spatial_coords(self):
    """
    Return spatial coordinates
    """
    metadata = self.metadata()
    x_c = {**dict(metadata.x().all_strings()), **dict(metadata.x().all_doubles())}
    y_c = {**dict(metadata.y().all_strings()), **dict(metadata.y().all_doubles())}

    s_c = {"x": x_c, "y": y_c}
    if metadata.n_spatial_dimensions() > 2:
        z_c = {**dict(metadata.z().all_strings()), **dict(metadata.z().all_doubles())}
        s_c.update({"z": z_c})
    return s_c

@property
def attrs(self):
    """
    Return attributes
    """
    metadata = self.metadata()    
    return {**dict(metadata.all_strings()), **dict(metadata.all_doubles())}

@property
def cf_mapping(self):
    """
    Return cf_mapping
    """
    mapping_info = self.grid().get_mapping_info()
    return {**dict(mapping_info.cf_mapping.all_strings()), **dict(mapping_info.cf_mapping.all_doubles())}
