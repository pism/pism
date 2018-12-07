def __getitem__(self, *args):
    return self.getitem(args[0][0], args[0][1])


def __setitem__(self, *args):
    if(len(args) == 2):
        self.setitem(args[0][0], args[0][1], args[1])
    else:
        raise ValueError("__setitem__ requires 2 arguments; received %d" % len(args))


def imshow(self, **kwargs):
    "Plot a 2D field using matplotlib.pylab.imshow()."

    try:
        import matplotlib.pylab as plt
    except:
        raise RuntimeError("Failed to import matplotlib.pylab. Please make sure that matplotlib is installed!")

    m = plt.imshow(self.numpy(), origin="lower", **kwargs)
    plt.colorbar(m)
    md = self.metadata()
    plt.title("{}, {}".format(md.get_string("long_name"), md.get_string("units")))
    plt.axis("equal")
