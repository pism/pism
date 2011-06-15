#!/usr/bin/env python

# Copyright (C) 2009-2010 Andy Aschwanden and Ed Bueler
# PISM SeaRISE Greenland worked example

#  this script fixes the mask in a PISM output to conform to mask choices
#    for the SeaRISE assessment
#  see:   http://websrv.cs.umt.edu/isis/index.php/Output_Format#Two-dimensional_output_variables
#  run postprocess.sh to use this one; see comments there
#  example usage (as done by postprocess.sh):
#    $ ./postprocess_mask.py UAF1_G_D3_C2_E0.nc

try:
    from netCDF4 import Dataset as CDF
except:
    from netCDF3 import Dataset as CDF
import sys

input = sys.argv[1]
print "  (postprocess_mask.py) fixing mask in %s to conform to SeaRISE ..." % input 

ice_thickness_threshold = 10.0

# PISM mask values:
pism_ice_free_bedrock = 0
pism_grounded = 2
pism_floating = 3
pism_ocean    = 4

# SeaRISE mask values:
searise_ocean = 0
searise_grounded_ice = 1
searise_floating_ice = 2
searise_ice_free_land = 3

nc = CDF(input, 'a')

# number of records to process:
try:
    N = len(nc.variables['t'][:])
except:
    N = len(nc.variables['time'][:])
    
thk = nc.variables['thk']
mask = nc.variables['mask']

mask.flag_meanings = "ocean grounded_ice floating_ice ice_free_land" ;
mask.flag_values = [0, 1, 2, 3];
mask.long_name = "integer mask specifying cell type"

for j in range(N):
    Mask = mask[j,:,:]
    Thk  = thk[j,:,:]
    tmp  = mask[j,:,:]

    # Combine PISM's grounded and dragging ice:
    tmp[(Mask == pism_grounded)] = searise_grounded_ice
    # Mark ice-free land:
    tmp[(Mask == pism_ice_free_bedrock)] = searise_ice_free_land

    # Mark ocean:
    tmp[Mask == pism_ocean] = searise_ocean
    # Mark floating ice:
    tmp[Mask == pism_floating] = searise_floating_ice

    mask[j,:,:] = tmp

# finish by updating history global attribute
import time
historysep = ' '
historystr = time.asctime() + ': ' + historysep.join(sys.argv) + '\n'
if 'history' in nc.ncattrs():
  nc.history = historystr + nc.history
else:
  nc.history = historystr

nc.close()

