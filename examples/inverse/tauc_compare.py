#!/usr/bin/env python3

try:
    import xarray as xr
except ImportError:
    print("xarray is not installed!")
    exit(1)
from matplotlib import pyplot as pp
from matplotlib import colors as mc
from optparse import OptionParser
from siple.reporting import endpause
import numpy as np

usage = """Usage: %prog [options]

Example: %prog -N 100 -n 0.1"""

parser = OptionParser(usage=usage)
parser.add_option("-i", "--input_file", type='string',
                  help='input file')
parser.add_option("-c", "--tauc_cap", type='float', default=200000,
                  help='maximum tauc value to display')

parser.add_option("-e", "--tauc_error_cap", type='float', default=0.2,
                  help='maximum relative error to display')


(options, args) = parser.parse_args()

try:
    ds = xr.open_dataset(options.input_file, decode_times=False, decode_cf=False)
except Exception:
    print('ERROR: option -i is required')
    parser.print_help()
    exit(0)

secpera = 3.15569259747e7
tauc = np.asarray(ds['tauc'].values).squeeze()
tauc_true = np.asarray(ds['tauc_true'].values).squeeze()

tauc_diff = tauc - tauc_true

not_ice = abs(np.asarray(ds['mask'].values).squeeze() - 2) > 0.01
tauc[not_ice] = 0
tauc_true[not_ice] = 0
tauc_diff[not_ice] = 0.


u_computed = np.asarray(ds['u_computed'].values).squeeze() * secpera
v_computed = np.asarray(ds['v_computed'].values).squeeze() * secpera
velbase_mag_computed = np.sqrt(u_computed * u_computed + v_computed * v_computed)

not_sliding = np.logical_and((abs(u_computed) < 10.), (abs(v_computed) < 10.))
tauc[not_ice] = 0
tauc_true[not_ice] = 0
tauc_diff[not_sliding] = 0.


# difference figure
pp.clf()
pp.imshow(tauc_diff.transpose() / tauc_true.transpose(), origin='lower',
          vmin=-options.tauc_error_cap, vmax=options.tauc_error_cap)
pp.title(r'$(\tau_c$ - true) / true')
pp.colorbar()

# side-by-side comparison
pp.figure()
pp.subplot(1, 2, 1)
pp.imshow(tauc.transpose(), origin='lower', vmin=0.0, vmax=options.tauc_cap)
pp.title(r'$\tau_c$  [from inversion]')
pp.colorbar()
pp.subplot(1, 2, 2)
pp.imshow(tauc_true.transpose(), origin='lower', vmin=0.0, vmax=options.tauc_cap)
pp.title(r'true $\tau_c$   [prior]')
pp.colorbar()

# show computed sliding speed
pp.figure()
im = pp.imshow(velbase_mag_computed.transpose(), origin='lower',
               norm=mc.LogNorm(vmin=0.1, vmax=1000.0))
pp.title('computed sliding speed')
t = [0.1, 1.0, 10.0, 100.0, 1000.0]
pp.colorbar(im, ticks=t, format='$%.1f$')

# pp.ion()
pp.show()
# endpause()
