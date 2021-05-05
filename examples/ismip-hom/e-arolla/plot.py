import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import netCDF4

# Plots PISM results for the ISMIP-HOM Experiment E.

f_sliding = netCDF4.Dataset("output-sliding.nc")
f_no_slip = netCDF4.Dataset("output-no-slip.nc")

# Surface ice velocity

x = f_sliding.variables["x"][:]
uvelsurf_sliding = f_sliding.variables["uvelsurf"][0,1,:]
uvelsurf_no_slip = f_no_slip.variables["uvelsurf"][0,1,:]

fig, ax = plt.subplots()
fig.dpi = 200
ax.plot(x, uvelsurf_no_slip, label="no slip")
ax.plot(x, uvelsurf_sliding, label="sliding")

ax.set(xlabel='x (meters)', ylabel='surface velocity (m/year) ',
       title='ISMIP-HOM Experiment E (Haut Glacier d’Arolla)')
ax.grid()
ax.legend()
fig.savefig("uvelsurf.png")

# Ice velocity throughout the ice volume
year = 365.2524 * 86400

topg = f_sliding.variables["topg"][0,1,:]
thk = f_sliding.variables["thk"][0,1,:]

uvel_sliding = f_sliding.variables["uvel_sigma"][0,1,:,:]
uvel_sliding *= year

Mx, Mz = uvel_sliding.shape

xx = np.repeat(x, Mz).reshape(Mx, Mz).T

zz = np.repeat(np.linspace(0, 1, Mz), Mx).reshape(Mz, Mx)
for k in range(Mz):
    zz[k, :] *= thk
    zz[k, :] += topg

# Sliding case
fig, ax = plt.subplots()
fig.dpi = 200

ax.set(xlabel='x (m)', ylabel='elevation (m)',
       title='ISMIP-HOM Exp E (sliding) ice velocity, m/year')

m = ax.pcolormesh(xx, zz, uvel_sliding.T, shading="nearest")
fig.colorbar(m)

fig.savefig("uvel_sliding.png")

# No-slip case
fig, ax = plt.subplots()
fig.dpi = 200

ax.set(xlabel='x (m)', ylabel='elevation (m)',
       title='ISMIP-HOM Exp E (no slip) ice velocity, m/year')

uvel_no_slip = f_no_slip.variables["uvel_sigma"][0,1,:,:]
uvel_no_slip *= year

m = ax.pcolormesh(xx, zz, uvel_no_slip.T, shading="nearest")
fig.colorbar(m)

fig.savefig("uvel_no_slip.png")
