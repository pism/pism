#!/usr/bin/env python3

from numpy import *
from matplotlib.pyplot import *
from sys import exit

try:
    import xarray as xr
except ImportError:
    print("xarray is not installed!")
    exit(1)

nameroot = "routing"

for dx in ("100", "50", "25", "15", "10", "5"):
    basename = nameroot + dx + "km"
    filename = basename + ".nc"
    print("%s: looking for file ..." % filename)
    try:
        nc = xr.open_dataset(filename, decode_times=False, decode_cf=False)
    except Exception:
        print("  can't read from file ...")
        continue

    x = asarray(squeeze(nc["x"].values))
    y = asarray(squeeze(nc["y"].values))

    for varname in ("bwat", "bwp", "psi"):   # psi must go after bwat, bwp
        print("  %s:  generating pcolor() image ..." % varname)
        lookup = "topg" if varname == "psi" else varname
        if lookup not in nc.variables:
            print("variable '%s' not found ... continuing ..." % varname)
            continue
        data = asarray(squeeze(nc[lookup].values)).transpose()

        if varname == "bwat":
            bwatdata = data.copy()
        if varname == "bwp":
            bwpdata = data.copy()

        if varname == "psi":
            # psi = bwp + rho_w g (topg + bwat)
            data = bwpdata + 1000.0 * 9.81 * (data + bwatdata)

        if varname == "bwat":
            units = "m"
            barmin = 0.0
            barmax = 650.0
            scale = 1.0
        else:
            units = "bar"
            barmin = -20.0
            barmax = 360.0
            scale = 1.0e5

        print("       [stats:  max = %9.3f %s, av = %8.3f %s]" %
              (data.max() / scale, units, data.sum() / (scale * x.size * y.size), units))
        pcolor(x / 1000.0, y / 1000.0, data / scale, vmin=barmin, vmax=barmax)
        colorbar()
        gca().set_aspect('equal')
        gca().autoscale(tight=True)
        xlabel('x  (km)')
        ylabel('y  (km)')
        dxpad = "%03d" % int(dx)
        pngfilename = varname + "_" + dxpad + "km" + ".png"
        print("    saving figure in %s ..." % pngfilename)
        savefig(pngfilename, dpi=300, bbox_inches='tight')
        close()

    nc.close()
