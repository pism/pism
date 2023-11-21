import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import sys
# Plots PISM results for the ISMIP-HOM Experiment E.

sliding = sys.argv[1]
no_slip = sys.argv[2]

f_sliding = netCDF4.Dataset(sliding)
f_no_slip = netCDF4.Dataset(no_slip)

# Surface ice velocity

x = f_sliding.variables["x"][:]
uvelsurf_sliding = f_sliding.variables["uvelsurf"][:]
uvelsurf_no_slip = f_no_slip.variables["uvelsurf"][:]

fig, ax = plt.subplots()
fig.dpi = 100
fig.set_size_inches(8, 4)

ax.plot(x, uvelsurf_no_slip, label="no slip")
ax.plot(x, uvelsurf_sliding, label="sliding")

ax.set(xlabel='x (meters)', ylabel='surface velocity (m/year) ')
ax.grid()
ax.legend()
fig.savefig("uvelsurf.png")

# Ice velocity throughout the ice volume
year = 365.2524 * 86400

topg = f_sliding.variables["topg"][:]
thk = f_sliding.variables["thk"][:]

def plot_uvel(data, thk, topg, title, filename):
    Mz, Mx = data.shape

    xx = np.repeat(x, Mz).reshape(Mx, Mz).T

    zz = np.repeat(np.linspace(0, 1, Mz), Mx).reshape(Mz, Mx)
    for k in range(Mz):
        zz[k, :] *= thk
        zz[k, :] += topg

    fig, ax = plt.subplots()

    # note: ignore the title (for now)
    ax.set(xlabel='x (m)', ylabel='elevation (m)')

    m = ax.pcolormesh(xx, zz, data, shading="gouraud")
    fig.colorbar(m, pad=0.15, orientation="horizontal", aspect=40)

    exaggeration = 4
    width = 8
    ratio = exaggeration * 700 / 5000
    fig.set_size_inches(width, ratio * width)
    fig.tight_layout()

    fig.dpi = 100
    fig.savefig(filename)

# Sliding case
uvel_sliding = f_sliding.variables["uvel_sigma"][:,:]
uvel_sliding *= year

plot_uvel(uvel_sliding, thk, topg,
          'ISMIP-HOM Exp E (sliding) ice velocity, m/year',
          "uvel_sliding.png")

# No-slip case
uvel_no_slip = f_no_slip.variables["uvel_sigma"][:,:]
uvel_no_slip *= year

plot_uvel(uvel_no_slip, thk, topg,
          'ISMIP-HOM Exp E (no slip) ice velocity, m/year',
          "uvel_no_slip.png")
