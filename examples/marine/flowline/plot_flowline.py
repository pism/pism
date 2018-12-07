#!/usr/bin/env python
# Copyright (C) 2011, 2014 Torsten Albrecht (torsten.albrecht@pik-potsdam.de)

import numpy as np
from pylab import figure, clf, hold
import matplotlib.pyplot as plt
try:
    from netCDF4 import Dataset as NC
except:
    print("netCDF4 is not installed!")
    sys.exit(1)

pism_output = "flowline_result_50km.nc"
pism_output = "flowline_result_20km.nc"
#pism_output = "flowline_result_10km.nc"
#pism_output = "flowline_result_5km.nc"
#pism_output = "flowline_result_2km.nc"
#pism_output = "flowline_result_1km.nc"

# load the PISM output
try:
    print("Loading PISM output from '%s'..." % (pism_output), end=' ')
    infile = NC(pism_output, 'r')

except Exception:
    print("""ERROR!\nSpecify NetCDF file from PISM run with -p.""")
    exit(-1)

printfile = pism_output[:-3] + ".pdf"

thk = np.squeeze(infile.variables["thk"][:])
ubar = np.squeeze(infile.variables["ubar"][:])

# "typical constant ice parameter" as defined in the paper and in Van der
# Veen's "Fundamentals of Glacier Dynamics", 1999
U0 = 300  # m/yr
H0 = 600  # m
spa = 3.1556926e7
# spa=3600*365.25*24
C = 2.4511e-18
#C = (rho_ice * standard_gravity * (1.0 - rho_ice/rho_ocean) / (4 * B0))**3
Hcf = 250  # calving thickness

#grid and time
x = np.squeeze(infile.variables["x"][:])
dx = x[1] - x[0]
t = np.squeeze(infile.variables["time"][:])
t = int(np.round(t / spa))

dist = 50.0
idist = int(dist / (dx / 1000.0))
H = thk[1, idist:]
u = ubar[1, idist:]
xd = x[idist:] / 1000.0 - dist


# exact van der veen solution
xcf = (Hcf ** -4 - H0 ** -4) * U0 * H0 / (4 * C * spa) * 0.001  # find 250m front position on 5km grid
ucf = U0 * H0 * (xcf * 1000.0 * 4 * C / (U0 * H0 / spa) + H0 ** (-4)) ** (0.25)  # exact (maximal) velocity at front

dxan = 1.0  # km
xan = np.arange(0, (x[len(x) - 1]) / 1000.0, dxan)
Han = np.zeros(len(xan))
Uan = np.zeros(len(xan))
for k in range(len(xan)):
    if xan[k] <= xcf:
        Han[k] = (k * dxan * 1000.0 * 4 * C / (U0 * H0 / spa) + H0 ** (-4)) ** (-0.25)
        Uan[k] = U0 * H0 * (k * dxan * 1000.0 * 4 * C / (U0 * H0 / spa) + H0 ** (-4)) ** (0.25)


# imshow plot
plt.rcParams['font.size'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['contour.negative_linestyle'] = 'solid'

plt.clf()
figure(1, figsize=(10, 9))
clf()
hold(True)
plt.subplots_adjust(left=0.14, bottom=0.14, right=0.86, top=None, wspace=0.08,
                    hspace=0.15)
plt.clf()
title = 'steady state ice thickness and velocity for ' + str(int(dx / 1000.0)) + ' km after ' + str(t) + ' yrs'
plt.suptitle(title)
# plot
ax1 = plt.subplot(211)
ax1.plot(xd, H, '#B26415', xan, Han, '#4E238B', linewidth=2)
leg1 = ax1.legend(('H_pism in m', 'H_exact in m'), 'upper right', shadow=True)
plt.xlim([0, 400])
# plt.xticks([])
# plt.xticks(np.arange(0,50,10))
plt.ylim([0, 600])
plt.grid(True)

arrowtext2a = 'max(u_exact)=%d m/yr at %d km' % (int(ucf), int(xcf))
arrowtext2b = 'max(u_pism)=%d m/yr at %d km' % (int(np.max(u)), int(xd[np.argmax(u)]))
ax2 = plt.subplot(212)
ax2.plot(xd, u, '#B26415', xan, Uan, '#4E238B', linewidth=2)
leg2 = ax2.legend((arrowtext2b, arrowtext2a), 'lower right', shadow=True)
plt.xlim([0, 400])
# plt.xticks([])
plt.ylim([0, 900])
plt.yticks(np.arange(0, 1000, 200))
# ax2.annotate(arrowtext2a, xy=(int(x[np.argmax(Uexact)]),int(np.max(Uexact))),  xycoords='data',
#             xytext=(40, 20), textcoords='offset points',
#             arrowprops=dict(arrowstyle="->"))
#
# ax2.annotate(arrowtext2b, xy=(int(x[np.argmax(u)]),int(np.max(u))),  xycoords='data',
#            xytext=(40, 0), textcoords='offset points', arrowprops=dict(arrowstyle="->"))
plt.grid(True)
plt.savefig(printfile, format='pdf')
plt.show()
