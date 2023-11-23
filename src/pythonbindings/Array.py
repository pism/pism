# code extending the Array class

def regrid(self, filename, critical=False, default_value=0.0):
    if critical == True:
        self._regrid(filename, Default.Nil())
    else:
        self._regrid(filename, Default(default_value))

def numpy(self):
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
