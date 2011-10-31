try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
from matplotlib import pyplot as pp
from matplotlib import colors as mc
from optparse import OptionParser
from siple.reporting import endpause
import PISM
import numpy as np

usage =  """Usage: %prog [options]

Example: %prog -N 100 -n 0.1"""

parser = OptionParser(usage=usage)
parser.add_option("-i","--input_file",type='string',
                    help='input file')
parser.add_option("-c","--cap",type='float',default=None,
                    help='maximum cap for difference')


(options, args) = parser.parse_args()

ds = netCDF.Dataset(options.input_file)

tauc = ds.variables['tauc'][...].squeeze()
tauc_true = ds.variables['tauc_true'][...].squeeze()

tauc_diff = tauc-tauc_true

not_ice = abs(ds.variables['mask'][...].squeeze() -2 ) > 0.01
tauc[not_ice] = 0
tauc_true[not_ice] = 0
tauc_diff[not_ice] = 0.


u_computed = ds.variables['u_computed'][...].squeeze()*PISM.secpera
v_computed = ds.variables['v_computed'][...].squeeze()*PISM.secpera
cbase_computed = np.sqrt(u_computed * u_computed + v_computed * v_computed)

not_sliding = np.logical_and( (abs(u_computed) < 10.) , (abs(v_computed) < 10.) )
tauc[not_ice]  = 0
tauc_true[not_ice] = 0
tauc_diff[not_sliding] = 0.

if not options.cap is None:
  C = options.cap
  tauc_diff = np.maximum(tauc_diff,-C)
  tauc_diff = np.minimum(tauc_diff,C)

  # tauc = np.maximum(tauc,-10*C)
  # tauc = np.minimum(tauc,10*C)
  # 
  # tauc_true = np.maximum(tauc_true,-10*C)
  # tauc_true = np.minimum(tauc_true,10*C)

# difference figure  
pp.clf()
pp.imshow(tauc_diff.transpose()/tauc_true.transpose(),origin='lower')
pp.title(r'$(\tau_c$ - true) / true')
pp.colorbar()

# side-by-side comparison like from 'vel2tauc.py -inv_final_plot'
pp.figure()
pp.subplot(1,2,1)
pp.imshow(tauc.transpose(),origin='lower',vmin=0.0,vmax=200000.0)
pp.title(r'$\tau_c$  [from inversion]')
pp.colorbar()
pp.subplot(1,2,2)
pp.imshow(tauc_true.transpose(),origin='lower',vmin=0.0,vmax=200000.0)
pp.title(r'true $\tau_c$   [prior]')
pp.colorbar()

# show computed sliding speed
pp.figure()
im = pp.imshow(cbase_computed.transpose(),origin='lower',
               norm=mc.LogNorm(vmin=0.1, vmax=1000.0))
pp.title('computed sliding speed')
t = [0.1, 1.0, 10.0, 100.0, 1000.0]
pp.colorbar(im, ticks=t, format='$%.1f$')

# pp.ion()
pp.show()
# endpause()
