import sys
from matplotlib import pyplot as plt
import numpy as np
import orbital_parameters as berger1978

input_file = sys.argv[1]
output_file = sys.argv[2]

# Compares orbital parameters from
#
# J. Laskar, P. Robutel, F. Joutel, M. Gastineau, A. C. M. Correia,
# and B. Levrard, “A long-term numerical solution for the insolation
# quantities of the Earth,” Astronomy &amp; Astrophysics, vol. 428,
# no. 1, pp. 261–285, Nov. 2004, doi: 10.1051/0004-6361:20041335.
#
# downloaded from https://vo.imcce.fr/insola/earth/online/earth/La2004/INSOLN.LA2004.BTL.ASC
#
# to ones computed using trigonometric expasions from
#
# A. L. Berger, “Long-term variations of daily insolation and
# quaternary climatic changes,” Journal of the Atmospheric Sciences,
# vol. 35, no. 12, pp. 2362–2367, Dec. 1978, doi:
# 10.1175/1520-0469(1978)035<2362:ltvodi>2.0.co;2.
#
# as implemented in orbital_parameters.py based on the code from NASA
# GISS ModelE.

# Laskar2004 version:
data = np.loadtxt(input_file)

# Laskar2004 uses units of "1000 years since J2000"
t                    = data[:, 0] * 1000 + 2000
eccentricity         = data[:, 1]
obliquity            = np.rad2deg(data[:, 2])
perihelion_longitude = np.rad2deg(data[:, 3])

fig, ax = plt.subplots(3, sharex=True)
fig.suptitle("Comparing orbital parameters from Laskar2004 and Berger1978")

# convert perihelion longitude from Laskar2004 from the heliocentric
# to the geocentric reference frame:
perihelion_longitude[:] += 180
perihelion_longitude[perihelion_longitude > 360] -= 360
perihelion_longitude[perihelion_longitude < 0]   += 360

ax[0].plot(t, eccentricity, label="Laskar2004")
ax[1].plot(t, obliquity, label="Laskar2004")
ax[2].plot(t, perihelion_longitude, label="Laskar2004 (converted to geocentric)")

# Berger 1978 version:

data2 = np.array([berger1978.orbital_parameters(y, restrict=True) for y in t])

eccentricity         = data2[:, 0]
obliquity            = np.rad2deg(data2[:, 1])
perihelion_longitude = np.rad2deg(data2[:, 2])

label="Berger1978 (orbital_parameters.py; geocentric)"
ax[0].plot(t, eccentricity, label=label)
ax[1].plot(t, obliquity, label=label)
ax[2].plot(t, perihelion_longitude, label=label)

x_min = -2.5e5

titles = ["eccentricity", "obliquity", "perihelion longitude"]
ylabel = ["1", "degree", "degree"]
for j in range(3):
    ax[j].set_xlim(xmin=x_min, xmax=2000)
    ax[j].grid(True)
    ax[j].set_title(titles[j])
    ax[j].legend()
    ax[j].set_ylabel(ylabel[j])

ax[2].set_xlabel("time relative to year 1, kyr")
fig.set_size_inches(20, 10)
fig.savefig(output_file, dpi=100)
