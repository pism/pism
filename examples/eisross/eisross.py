#! /usr/bin/env python

import sys
import getopt
import time
from numpy import *
try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

# set constants
SECPERA = 3.1556926e7
GRID_FILE = '111by147Grid.dat'
KBC_FILE = 'kbc.dat'
INLETS_FILE = 'inlets.dat'
WRIT_FILE = 'ross.nc'
dxROSS = 6822.0 # meters
MASK_BC    = -2
MASK_SHEET =  1
MASK_FLOATING = 3
VERBOSE = 0

# function which will print ignored lines if VERBOSE > 0
def vprint(s):
    if VERBOSE > 0:
        print s

# function to read a 2d variable from EISMINT-ROSS data file
# allows choice of _FillValue, shifting, and scaling
def read2dROSSfloat(mygrid,myarray,xs,xm,My,mymissing,myshift,myscale):
    vprint(mygrid.readline()) # ignore two lines
    vprint(mygrid.readline())
    Mx = xs + xm
    for i in range(xs):
        for j in range(My):
            myarray[i,j] = mymissing
    for i in range(xm):
        j = 0
        for num in mygrid.readline().split():
            myarray[i+xs,j] = (float(num) + myshift) * myscale
            j = j + 1

# function convert velocities from (azimuth,magnitude), with magnitude to (u,v)
def uvGet(mag,azi):
    u = mag * sin((pi/180.0) * azi)
    v = mag * cos((pi/180.0) * azi)
    return (u,v)

##### command line arguments #####
try:
    opts, args = getopt.getopt(sys.argv[1:], "p:o:v:", ["prefix=", "out=", "verbose="])
    for opt, arg in opts:
        if opt in ("-p", "--prefix"):
            GRID_FILE = arg + GRID_FILE
            KBC_FILE = arg + KBC_FILE
            INLETS_FILE = arg + INLETS_FILE
        if opt in ("-o", "--out"):
            WRIT_FILE = arg
        if opt in ("-v", "--verbose"):
            verbose = float(arg)
except getopt.GetoptError:
    print 'Incorrect command line arguments'
    sys.exit(2)


##### read 111by147Grid.dat #####
print "reading grid data from ",GRID_FILE
vprint("VERBOSE > 0, so printing ignored lines in " + GRID_FILE)
grid=open(GRID_FILE, 'r')
vprint(grid.readline()) # ignore first line
# second line gives dimensions; read and allocate accordingly
dim=[]
for num in grid.readline().split():
    dim.append(num)
xmROSS = int(dim[0])
MyROSS = int(dim[1])
xsROSS = MyROSS - xmROSS
MxROSS = MyROSS
# compute RIGGS grid coordinates (from original RIGGS documents)
lat = zeros((MxROSS, MyROSS), float32)
lon = zeros((MxROSS, MyROSS), float32)
dlat = (-5.42445 - (-12.3325)) / 110.0
dlon = (3.72207 - (-5.26168)) / 146.0
for i in range(MxROSS):
    for j in range(MyROSS):
        lat[i,j] = -12.3325 - dlat * 46.0 + i * dlat
        lon[i,j] = -5.26168 + j * dlon
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
vprint(grid.readline()) # ignore two more lines
vprint(grid.readline())
j=0;
for line in range(xsROSS):
    for i in range(MyROSS):
        eislat[j,i] = 9999.;
    j = j + 1
for line in range(xmROSS):
    latvalue = float(grid.readline())
    for i in range(MyROSS):
        eislat[j,i] = latvalue
    j = j + 1
vprint(grid.readline()) # read extra value
# note there are actually 148 "columns position" values in 111by147Grid.dat file
vprint(grid.readline()) # ignore two lines
vprint(grid.readline())
i=0;
for line in range(MyROSS):
    lonvalue = float(grid.readline())
    for j in range(MxROSS):
        eislon[j,i] = lonvalue
    i = i + 1
vprint(grid.readline()) # read extra value
vprint(grid.readline()) # ignore two lines
vprint(grid.readline())
for i in [0, 1]:
    for j in range(MyROSS):
        mask[i,j] = MASK_BC
for i in range(xsROSS-2):
    for j in range(MyROSS):
        mask[i+2,j] = MASK_FLOATING
for i in range(xmROSS):
    j = 0
    for num in grid.readline().split():
        if int(num) == 1:
            mask[i+xsROSS,j] = MASK_FLOATING
        else:
            mask[i+xsROSS,j] = MASK_BC
        j = j + 1
read2dROSSfloat(grid,azi,xsROSS,xmROSS,MyROSS,9999.,0.0,1.0)
read2dROSSfloat(grid,mag,xsROSS,xmROSS,MyROSS,9999.,0.0,1.0 / SECPERA)
read2dROSSfloat(grid,thk,xsROSS,xmROSS,MyROSS,1.0,0.0,1.0)
vprint(grid.readline()) # ignore two lines
vprint(grid.readline())
for i in range(xsROSS):
    for j in range(MyROSS):
        accur[i,j] = -1
for i in range(xmROSS):
    j = 0
    for num in grid.readline().split():
        if float(num) == 1.0:
            accur[i+xsROSS,j] = 1
        else:
            accur[i+xsROSS,j] = 0
        j = j + 1
read2dROSSfloat(grid,bed,xsROSS,xmROSS,MyROSS,-600.0,0.0,-1.0)
# set thickness to 1.0 m according to this info
vprint(grid.readline()) # ignore two lines
vprint(grid.readline())
for i in range(xmROSS):
    j = 0
    for num in grid.readline().split():
        if (float(num) == 1.0):
            thk[i,j] = 1.0
        j = j + 1
read2dROSSfloat(grid,accum,xsROSS,xmROSS,MyROSS,0.2/SECPERA,0.0,1.0/(SECPERA * 1000.0))
read2dROSSfloat(grid,barB,xsROSS,xmROSS,MyROSS,9999.,0.0,1.0)
read2dROSSfloat(grid,Ts,xsROSS,xmROSS,MyROSS,248.0,273.15,1.0)
grid.close()


##### create arrays for observed ubar, vbar and fill with _FillValue #####
ubarOBS = zeros((MxROSS, MyROSS), float32)
vbarOBS = zeros((MxROSS, MyROSS), float32)

##### read kbc.dat #####
print "reading boundary condition locations from ",KBC_FILE
kbc=open(KBC_FILE, 'r')
for count in range(77):
    coords = kbc.readline().split()
    i = int(coords[0]) + xsROSS
    j = int(coords[1])
    mask[i,j] = MASK_BC
    [ubarOBS[i,j], vbarOBS[i,j]] = uvGet(mag[i,j],azi[i,j])
kbc.close()

##### read inlets.dat #####
print "reading additional boundary condition locations and data"
print "   from ",INLETS_FILE
inlets=open(INLETS_FILE, 'r')
for count in range(22):
    data = inlets.readline().split()
    i = int(data[0]) + xsROSS
    j = int(data[1])
    mask[i,j] = MASK_BC
    [ubarOBS[i,j], vbarOBS[i,j]]  = uvGet(float(data[3]) / SECPERA,float(data[2]))
inlets.close()

##### compute coordinates #####
x = zeros(MxROSS)
y = zeros(MyROSS)
for i in range(MxROSS):
    x[i] = dxROSS * float(i - (MxROSS - 1)/2)
for j in range(MyROSS):
    y[j] = dxROSS * float(j - (MyROSS - 1)/2)

##### define dimensions in NetCDF file #####
ncfile = NC(WRIT_FILE, 'w',format='NETCDF3_CLASSIC')
xdim = ncfile.createDimension('y', MyROSS)
ydim = ncfile.createDimension('x', MxROSS)

##### define variables, set attributes, write data #####
# format: ['units', 'long_name', 'standard_name', '_FillValue', array]
vars = {'y': ['m',
              'x-coordinate in Cartesian system',
              'projection_x_coordinate',
              None,
              x],
        'x': ['m',
              'y-coordinate in Cartesian system',
              'projection_y_coordinate',
              None,
              y],
        'lat': ['degrees_north',
                'RIGGS grid south latitude',
                'latitude',
                None,
                lat],
        'lon': ['degrees_east',
                'RIGGS grid west longitude',
                'longitude',
                None,
                lon],
        'mask': [None,
                 'grounded or floating integer mask',
                 None,
                 None,
                 mask],
        'azi_obs': ['degrees_east',
                    'EISMINT ROSS observed ice velocity azimuth',
                    None,
                    9999.0,
                    azi],
        'mag_obs': ['m s-1',
                    'EISMINT ROSS observed ice velocity magnitude',
                    None,
                    9999.0,
                    mag],
        'thk': ['m',
                'floating ice shelf thickness',
                'land_ice_thickness',
                1.0,
                thk],
        'accur': [None,
                  'EISMINT ROSS flag for accurate observed velocity',
                  None,
                  -1,
                  accur],
        'topg': ['m',
                 'bedrock surface elevation',
                 'bedrock_altitude',
                 -600.0,
                 bed],
        'acab': ['m s-1',
                 'mean annual net ice equivalent accumulation rate',
                 'land_ice_surface_specific_mass_balance',
                 0.2/SECPERA,
                 accum],
        'barB': ['Pa^(1/3)',
                 'vertically-averaged ice hardness coefficient',
                 None,
                 9999.0,
                 barB],
        'artm': ['K',
                 'annual mean air temperature at ice surface',
                 'surface_temperature',
                 248.0,
                 Ts],
        'ubar': ['m s-1',
                 'vertical average of horizontal velocity of ice in projection_x_coordinate direction',
                 'land_ice_vertical_mean_x_velocity',
                 1/SECPERA,
                 ubarOBS],
        'vbar': ['m s-1',
                 'vertical average of horizontal velocity of ice in projection_y_coordinate direction',
                 'land_ice_vertical_mean_y_velocity',
                 1/SECPERA,
                 vbarOBS],}

for name in vars.keys():
    [_, _, _, fill_value, data] = vars[name]
    if name in ['x', 'y']:
        var = ncfile.createVariable(name, 'f4', (name,))
    else:
        var = ncfile.createVariable(name, 'f4', ('y', 'x'), fill_value = fill_value)

    for each in zip(['units', 'long_name', 'standard_name'], vars[name]):
        if each[1]:
            setattr(var, each[0], each[1])

    var[:] = data

##### attributes in NetCDF file #####
# set global attributes
historysep = ' '
historystr = time.asctime() + ': ' + historysep.join(sys.argv) + '\n'
setattr(ncfile, 'history', historystr)
setattr(ncfile, 'Conventions', 'CF-1.4') # only global attribute

# finish up
ncfile.close()
print "NetCDF file ",WRIT_FILE," created"

