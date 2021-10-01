#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import netCDF4 as nc

# ICE FRONT 

data_path = ""

experiments = [data_path+'result_1a_A1_800_300_1750e3_601_onesided_',
               data_path+'result_1a_A1_800_300_1750e3_601_centered_']
labels = ['one sided', 'centered']

# FIXME resolution dependent!!

refdat = data_path+'ICESHELF_1a_A1_800_300_1750e3_601_default.nc'
with nc.Dataset(refdat, 'r') as ncr:
    velref = np.ma.array(ncr.variables["u_ssa_bc"][1,:])

with nc.Dataset(experiments[0]+'default.nc', 'r') as ncr:
    mask = np.ma.array(ncr.variables["mask"][0,1,:])
    x = np.ma.array(ncr.variables["x"][:]) / 1000.0 # convert to km

start_idx = 2

# find the index of the last grounded cell gl_idx:
MASK_FLOATING = 3
gl_idx = 0
for i in range(len(x) - 1):
    if mask[i] < MASK_FLOATING and mask[i + 1] >= MASK_FLOATING:
        gl_idx = i
        break

#MISMIP
secpera=3.15569259747e7

def plot_taud(axis, exp, i):
    varname = "taud_x"
    axis.axvline(x=(x[gl_idx]), color='grey', linestyle='--', label=None)
    axis.axhline(y=0, color='grey', linestyle='--', label=None)
    for n,run in enumerate(experiments):

        with nc.Dataset(run+"default.nc", 'r') as ncr:
            refvariable = np.ma.array(ncr.variables[varname][0,1,:])

        with nc.Dataset(run+exp, 'r') as ncr:
            variable = np.ma.array(ncr.variables[varname][0,1,:])

        variable.mask = mask > MASK_FLOATING
        axis.plot(x[start_idx:], variable[start_idx:]-refvariable[start_idx:],
                  linewidth=2*(2-n),
                  label=labels[n])
        axis.scatter(x[start_idx:], variable[start_idx:]-refvariable[start_idx:], label=None)

    axis.set_ylim([-1200,600]) # km
    axis.set_xlim([1700,1755]) # km
    axis.text(1702,400,exp.split('_')[-1].rstrip('.nc'))

    if i==0:
      axis.set_title('Change in driving stress (Pa)')
    if i<3:
      axis.set_xticklabels([])

def plot_velocity(axis, exp, i):
    varname = "u_ssa"
    axis.axvline(x=x[gl_idx], color='grey', linestyle='--', label=None)
    axis.axhline(y=0, color='grey', linestyle='--', label=None)
    for n,run in enumerate(experiments):
        refvariable=velref

        with nc.Dataset(run+exp, 'r') as ncr:
            variable = np.ma.array(ncr.variables[varname][0,1,:])*secpera

        variable.mask = mask > MASK_FLOATING
        axis.plot(x[start_idx:], variable[start_idx:] - refvariable[start_idx:],
                  linewidth=2*(2-n),
                  label=labels[n])

    axis.set_ylim([-2000,500]) # km
    axis.set_xlim([0,1760]) # km

    if i==0:
        axis.set_title('Change in velocity (m/a)')
    if i<3:
        axis.set_xticklabels([])

fig, axes = plt.subplots(4, 2, sharex=False, sharey=False, figsize=(10, 8))
plt.subplots_adjust( wspace=0.25, hspace=0.1)

for i,exp in enumerate([ 'default.nc', 'p1.nc', 'p2.nc', 'p3.nc']):
    # plot changes is taud
    axis = axes.flatten()[2*i]
    plot_taud(axis, exp, i)

    # plot changes in ice velocity
    axis = axes.flatten()[2*i+1]
    plot_velocity(axis, exp, i)

ax = axes.flatten()[0]
ax.legend(loc='lower right')
axes.flatten()[6].set_xlabel('Distance from ice divide (km)')
axes.flatten()[7].set_xlabel('Distance from ice divide (km)')

plt.show()
