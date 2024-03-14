#!/usr/bin/env python3


# Creates circular setup fpr CalvingMIP as in https://github.com/JRowanJordan/CalvingMIP/wiki/Circular-domain
# Model parameters, see https://github.com/JRowanJordan/CalvingMIP/wiki/Experimental-parameters

# run as "python setup_circular.py -L 1.6e6 -M 321" for 5km resolution and 1600 x 1600 km domain width

# set path to ../../preprocessing/PISMNC.py
from PISMNC import PISMDataset as NC
from optparse import OptionParser
import numpy as np


R=800e3
Bc=900
Bl=-2000
Ba=1100
rc=0

cd=700e3
rhoi=910.0 #917.0
ac=0.3*rhoi
t0=-10.0
h0=100.0

parser = OptionParser()

parser.usage = "%prog [options]"
parser.description = "Fills missing values in variables selected using -v in the file given by -f."
parser.add_option("-o", dest="output_file_name", default="circular_input.nc",
                  help="output file name")
parser.add_option("-L", dest="length", default=100e3, type=float,
                  help="domain length, meters")
parser.add_option("-M", dest="Mx", default=101, type=int,
                  help="number of grid points in the flow direction")

(options, args) = parser.parse_args()


nc = NC(options.output_file_name, 'w')

x = np.linspace(-options.length/2, options.length/2, options.Mx)
xx, yy = np.meshgrid(x, x)
zeros = np.zeros((options.Mx, options.Mx))


r=np.sqrt(xx**2+yy**2)
theta=np.arctan2(yy,xx)
topg = Bc-(Bc-Bl)*(r-rc)**2./(R-rc)**2


zeros = np.zeros((options.Mx, options.Mx))
thk = zeros.copy() + h0
thk[r>=cd]=0.0

v = zeros.copy()
u = zeros.copy()

temp = np.ones_like(zeros)*t0   # C
smb  = np.ones_like(zeros)*ac   # m/yr*kg/m3

calvmask = np.ones_like(zeros)
calvmask[r>=cd]=0.0

nc.create_dimensions(x, x)
nc.write("topg", topg, attrs={"units": "m", "long_name": "bed_topography"})
nc.write("climatic_mass_balance", smb, attrs={"units": "kg m-2 year-1"})
nc.write("ice_surface_temp", temp, attrs={"units": "Celsius"})
nc.write("thk", thk,
         attrs={"units": "m", "standard_name": "land_ice_thickness"})
nc.write("land_ice_area_fraction_retreat", calvmask,
         attrs={"long_name": "maximum ice extent mask"})
nc.write("ubar", u,
         attrs={"units": "m/year", "long_name": "x-component of velocity"})
nc.write("vbar", v,
         attrs={"units": "m/year", "long_name": "y-component of velocity"})

nc.close()
