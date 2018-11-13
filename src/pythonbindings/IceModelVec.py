# code extending the IceModelVec class


def regrid(self, filename, critical=False, default_value=0.0):
    if critical == True:
        flag = CRITICAL
    else:
        flag = OPTIONAL
    self._regrid(filename, flag, default_value)


def numpy(self):
    "Return a NumPy array containing data from this field (on rank 0)."
    tmp = self.allocate_proc0_copy()
    self.put_on_proc0(tmp.get())
    import numpy

    if self.grid().ctx().rank() == 0:
        return numpy.array(tmp.get()).reshape(self.shape())
    else:
        return None
