def __getitem__(self, *args):
    return self.getitem(args[0][0], args[0][1])


def __setitem__(self, *args):
    if(len(args) == 2):
        i = args[0][0]
        j = args[0][1]
        val = args[1]
        if(isinstance(val, list) and len(val) == 2):
            self.setitem(i, j, val[0], val[1])
        else:
            self.setitem(i, j, val)
    else:
        raise ValueError("__setitem__ requires 2 arguments; received %d" % len(args))
