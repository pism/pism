def get(self, key):
    methods = [self.get_2d_cell_type,
               self.get_2d_scalar,
               self.get_2d_vector,
               self.get_3d_scalar]

    for method in methods:
        try:
            return method(key)
        except:
            pass

    return self._get(key)
