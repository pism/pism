#! /usr/bin/env python

import sys
import getopt
from numpy import *
from pycdf import *

# set constants
SECPERA = 3.1556926e7
GRID_FILE = '111by147Grid.dat'
KBC_FILE = 'kbc.dat'
INLETS_FILE = 'inlets.dat'
WRIT_FILE = 'ross.nc'
dxROSS = 6822.0 # meters
MASK_SHEET = 1
MASK_FLOATING = 3
VERBOSE = 0

# function which will print ignored lines if VERBOSE > 0
def vprint(s):
  if VERBOSE > 0:
    print s

# function to read a 2d variable from EISMINT-ROSS data file
# allows choice of missing_value, shifting, and scaling
def read2dROSSfloat(mygrid,myarray,xs,xm,My,mymissing,myshift,myscale):
  vprint(mygrid.readline()) # ignor two lines
  vprint(mygrid.readline())
  Mx = xs + xm
  top = Mx - 1
  for i in range(xs):
    for j in range(My):
      myarray[top - i,j] = mymissing
  for i in range(xm):
    j = 0
    for num in mygrid.readline().split():
      myarray[top - (i+xs),j] = (float(num) + myshift) * myscale
      j = j + 1

# function convert velocities from (azimuth,magnitude), with magnitude to (u,v)
def uvGet(mag,azi):
  uv = zeros((2,),float32)
  uv[0] = mag * cos((pi/180.0) * azi)
  uv[1] = mag * sin((pi/180.0) * azi)
  return uv

##### command line arguments #####
try:
  opts, args = getopt.getopt(sys.argv[1:], "p:o:", ["prefix=", "out="])
  for opt, arg in opts:
    if opt in ("-p", "--prefix"):
      GRID_FILE = arg + GRID_FILE
      KBC_FILE = arg + KBC_FILE
      INLETS_FILE = arg + INLETS_FILE
    if opt in ("-o", "--out"):
      WRIT_FILE = arg
except getopt.GetoptError:
  print 'Incorrect command line arguments'
  sys.exit(2)
  

##### read 111by147Grid.dat #####
print "reading grid data from ",GRID_FILE
vprint("VERBOSE > 0, so printing ignored lines in " + GRID_FILE)
grid=open(GRID_FILE, 'r')
vprint(grid.readline()) # ignor first line
# second line gives dimensions; read and allocate accordingly
dim=[]
for num in grid.readline().split():
   dim.append(num)
xmROSS = int(dim[0])
MyROSS = int(dim[1])
xsROSS = MyROSS - xmROSS
MxROSS = MyROSS
top = MxROSS - 1
# compute RIGGS grid coordinates (from original RIGGS documents)
lat = zeros((MxROSS, MyROSS), float32)
lon = zeros((MxROSS, MyROSS), float32)
dlat = (-5.42445 - (-12.3325)) / 110.0
dlon = (3.72207 - (-5.26168)) / 146.0
#from Matlab file ross_plot.m:
#dlat = (-5.42445 - (-12.3325))/110;
#gridlatext = linspace(-12.3325 - dlat * 46,-5.42445,147);
#gridlon = linspace(-5.26168,3.72207,147);
for i in range(MxROSS):
  for j in range(MyROSS):
    lat[top - i,j] = -12.3325 - dlat * 46.0 + i * dlat
    lon[top - i,j] = -5.26168 + j * dlon
#these are to be filled from 111by147.dat:
eislat = zeros((MxROSS, MyROSS), float32) # actually ignored
eislon = zeros((MxROSS, MyROSS), float32) # actually ignored
mask = zeros((MxROSS,MyROSS), int16)
azi = zeros((MxROSS, MyROSS), float32)
mag = zeros((MxROSS, MyROSS), float32)
thk = zeros((MxROSS, MyROSS), float32)
accur = zeros((MxROSS,MyROSS), int16)
bed = zeros((MxROSS, MyROSS), float32)
accum = zeros((MxROSS, MyROSS), float32)
barB = zeros((MxROSS, MyROSS), float32)
Ts = zeros((MxROSS, MyROSS), float32)
# note there are actually 112 "rows position" values in 111by147Grid.dat file
vprint(grid.readline()) # ignor two more lines
vprint(grid.readline())
j=0;
for line in range(xsROSS):
   for i in range(MyROSS):
      eislat[top-j,i] = 9999.;
   j = j + 1
for line in range(xmROSS):
   latvalue = float(grid.readline())
   for i in range(MyROSS):
      eislat[top-j,i] = latvalue
   j = j + 1
vprint(grid.readline()) # read extra value
# note there are actually 148 "columns position" values in 111by147Grid.dat file
vprint(grid.readline()) # ignor two lines
vprint(grid.readline())
i=0;
for line in range(MyROSS):
   lonvalue = float(grid.readline())
   for j in range(MxROSS):
      eislon[top-j,i] = lonvalue
   i = i + 1
vprint(grid.readline()) # read extra value
vprint(grid.readline()) # ignor two lines
vprint(grid.readline())
for i in [0, 1]:
  for j in range(MyROSS):
    mask[top - i,j] = MASK_SHEET
for i in range(xsROSS-2):
  for j in range(MyROSS):
    mask[top - (i+2),j] = MASK_FLOATING
for i in range(xmROSS):
  j = 0
  for num in grid.readline().split():
    if int(num) == 1:
       mask[top - (i+xsROSS),j] = MASK_FLOATING
    else:
       mask[top - (i+xsROSS),j] = MASK_SHEET
    j = j + 1
read2dROSSfloat(grid,azi,xsROSS,xmROSS,MyROSS,9999.,0.0,1.0)
read2dROSSfloat(grid,mag,xsROSS,xmROSS,MyROSS,9999.,0.0,1.0 / SECPERA)
read2dROSSfloat(grid,thk,xsROSS,xmROSS,MyROSS,1.0,0.0,1.0)
vprint(grid.readline()) # ignor two lines
vprint(grid.readline())
for i in range(xsROSS):
  for j in range(MyROSS):
    accur[top - i,j] = -1
for i in range(xmROSS):
  j = 0
  for num in grid.readline().split():
    if float(num) == 1.0:
       accur[top - (i+xsROSS),j] = 1
    else:
       accur[top - (i+xsROSS),j] = 0
    j = j + 1
read2dROSSfloat(grid,bed,xsROSS,xmROSS,MyROSS,-600.0,0.0,-1.0)
# set thickness to 1.0 m according to this info
vprint(grid.readline()) # ignor two lines
vprint(grid.readline())
for i in range(xmROSS):
  j = 0
  for num in grid.readline().split():
    if (float(num) == 1.0):
      thk[top - i,j] = 1.0
    j = j + 1
read2dROSSfloat(grid,accum,xsROSS,xmROSS,MyROSS,0.2/SECPERA,0.0,1.0/(SECPERA * 1000.0))
read2dROSSfloat(grid,barB,xsROSS,xmROSS,MyROSS,9999.,0.0,1.0)
read2dROSSfloat(grid,Ts,xsROSS,xmROSS,MyROSS,248.0,273.15,1.0)
grid.close()


##### create arrays for observed ubar, vbar and fill with missing_value #####
ubarOBS = zeros((MxROSS, MyROSS), float32)
vbarOBS = zeros((MxROSS, MyROSS), float32)
bcflag = zeros((MxROSS, MyROSS), int16)
for i in range(MxROSS):
  for j in range(MxROSS):
    ubarOBS[i,j] = 1.0 / SECPERA
    vbarOBS[i,j] = 1.0 / SECPERA
    bcflag[i,j] = 0
# also fill in zeros along sides; better for Laplace solution
for i in range(MxROSS):
  ubarOBS[i,0] = 0.0
  vbarOBS[i,0] = 0.0
  ubarOBS[i,MyROSS-1] = 0.0
  vbarOBS[i,MyROSS-1] = 0.0
for j in range(MyROSS):
  ubarOBS[0,j] = 0.0
  vbarOBS[0,j] = 0.0
  ubarOBS[MxROSS-1,j] = 0.0
  vbarOBS[MxROSS-1,j] = 0.0
  
##### read kbc.dat #####
print "reading boundary condition locations from ",KBC_FILE
kbc=open(KBC_FILE, 'r')
for count in range(77):
  coords = kbc.readline().split()
  i = top - (int(coords[0]) + xsROSS)
  j = int(coords[1])
  mask[i,j] = MASK_SHEET
  bcflag[i,j] = 1
  [ubarOBS[i,j], vbarOBS[i,j]] = uvGet(mag[i,j],azi[i,j])
kbc.close()

##### read inlets.dat #####
print "reading additional boundary condition locations and data"
print "   from ",INLETS_FILE
inlets=open(INLETS_FILE, 'r')
for count in range(22):
  data = inlets.readline().split()
  i = top - (int(data[0]) + xsROSS)
  j = int(data[1])
  mask[i,j] = MASK_SHEET
  bcflag[i,j] = 1
  [ubarOBS[i,j], vbarOBS[i,j]]  = uvGet(float(data[3]) / SECPERA,float(data[2]))
inlets.close()


##### create and define dimensions and variables in NetCDF file #####
ncfile = CDF(WRIT_FILE, NC.WRITE|NC.CREATE|NC.TRUNC)
ncfile.automode()
# define the dimensions
xdim = ncfile.def_dim('x', MxROSS)
ydim = ncfile.def_dim('y', MyROSS)
# define the variables
xvar = ncfile.def_var('x', NC.FLOAT, (xdim,))
yvar = ncfile.def_var('y', NC.FLOAT, (ydim,))
latvar = ncfile.def_var('lat', NC.FLOAT, (xdim,ydim))
lonvar = ncfile.def_var('lon', NC.FLOAT, (xdim,ydim))
maskvar = ncfile.def_var('mask', NC.INT, (xdim,ydim))
azivar = ncfile.def_var('azi_obs', NC.FLOAT, (xdim,ydim))
magvar = ncfile.def_var('mag_obs', NC.FLOAT, (xdim,ydim))
thkvar = ncfile.def_var('thk', NC.FLOAT, (xdim,ydim))
accurvar = ncfile.def_var('accur', NC.INT, (xdim,ydim))
bedvar = ncfile.def_var('topg', NC.FLOAT, (xdim,ydim))
accumvar = ncfile.def_var('acab', NC.FLOAT, (xdim,ydim))
barBvar = ncfile.def_var('barB', NC.FLOAT, (xdim,ydim))
Tsvar = ncfile.def_var('artm', NC.FLOAT, (xdim,ydim))
ubarvar = ncfile.def_var('ubar', NC.FLOAT, (xdim,ydim))
vbarvar = ncfile.def_var('vbar', NC.FLOAT, (xdim,ydim))
bcflagvar = ncfile.def_var('bcflag', NC.INT, (xdim,ydim))

##### attributes in NetCDF file #####
# set global attributes
setattr(ncfile, 'Conventions', 'CF-1.0')
# set the attributes of the variables
setattr(xvar, 'axis', 'X')
setattr(xvar, 'long_name', 'x-coordinate in Cartesian system')
setattr(xvar, 'standard_name', 'projection_x_coordinate')
setattr(xvar, 'units', 'm')
setattr(yvar, 'axis', 'Y')
setattr(yvar, 'long_name', 'y-coordinate in Cartesian system')
setattr(yvar, 'standard_name', 'projection_y_coordinate')
setattr(yvar, 'units', 'm')
setattr(latvar, 'long_name', 'RIGGS grid south latitude')
setattr(lonvar, 'long_name', 'RIGGS grid west longitude')
setattr(maskvar, 'long_name', 'grounded or floating integer mask')
setattr(azivar, 'long_name', 'EISMINT ROSS observed ice velocity azimuth')
setattr(azivar, 'units', 'degrees east')
setattr(azivar, 'missing_value', 9999.)
setattr(magvar, 'long_name', 'EISMINT ROSS observed ice velocity magnitude')
setattr(magvar, 'units', 'm s-1')
setattr(magvar, 'missing_value', 9999.)
setattr(thkvar, 'long_name', 'floating ice thickness')
setattr(thkvar, 'units', 'm')
setattr(thkvar, 'missing_value', 1.)
setattr(accurvar, 'long_name', 'EISMINT ROSS flag for accurate observed velocity')
setattr(accurvar, 'missing_value', -1)
setattr(bedvar, 'long_name', 'bedrock surface elevation')
setattr(bedvar, 'standard_name', 'bedrock_altitude')
setattr(bedvar, 'units', 'm')
setattr(bedvar, 'missing_value', -600.0)
setattr(accumvar, 'long_name', 'mean annual net ice equivalent accumulation rate')
setattr(accumvar, 'standard_name', 'land_ice_surface_specific_mass_balance')
setattr(accumvar, 'units', 'm s-1')
setattr(accumvar, 'missing_value', 0.2 / SECPERA)
setattr(barBvar, 'long_name', 'vertically-averaged ice hardness coefficient')
setattr(barBvar, 'units', 'Pa^(1/3)')
setattr(barBvar, 'missing_value', 9999.)
setattr(Tsvar, 'long_name', 'annual mean air temperature at ice surface')
setattr(Tsvar, 'standard_name', 'surface_temperature')
setattr(Tsvar, 'units', 'K')
setattr(Tsvar, 'missing_value', 248.0)
setattr(ubarvar, 'long_name', 
        'vertical average of horizontal velocity of ice in projection_x_coordinate direction')
setattr(ubarvar, 'units', 'm s-1')
setattr(ubarvar, 'missing_value', 1.0 / SECPERA)
setattr(vbarvar, 'long_name', 
        'vertical average of horizontal velocity of ice in projection_y_coordinate direction')
setattr(vbarvar, 'units', 'm s-1')
setattr(vbarvar, 'missing_value', 1.0 / SECPERA)
setattr(bcflagvar, 'long_name', 'location of Dirichlet boundary condition for velocity')

##### write data into and close NetCDF file #####
for i in range(MxROSS):
	xvar[i] = dxROSS * float(i - (MxROSS - 1)/2)
for j in range(MyROSS):
	yvar[j] = dxROSS * float(j - (MyROSS - 1)/2)
latvar[:] = lat
lonvar[:] = lon
maskvar[:] = mask
azivar[:] = azi
magvar[:] = mag
thkvar[:] = thk
accurvar[:] = accur
bedvar[:] = bed
accumvar[:] = accum
barBvar[:] = barB
Tsvar[:] = Ts
ubarvar[:] = ubarOBS
vbarvar[:] = vbarOBS
bcflagvar[:] = bcflag
# finish up
ncfile.close()
print "NetCDF file ",WRIT_FILE," created"

