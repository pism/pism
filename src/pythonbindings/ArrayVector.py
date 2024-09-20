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

