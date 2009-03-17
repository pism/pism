#!/usr/bin/env python

# solverSM.py is a Matlab-to-pylab translation of "SMsolver2.m" and other
# Matlab codes contained in
#    http://homepages.ulb.ac.be/~fpattyn/mismip/MISMIP_distribution.tar
# see major comments at the start of SMsolver2.m; comments here 
# starting "#SM" are from those Matlab codes

from scipy import *
from pylab import *
import sys
import time
from getopt import getopt, GetoptError

# specify which MISMIP using options, e.g.
#   ./solverSM.py --exper=3 --sliding=a --step=8 --out=foo.png
# or equivalently
#   ./solverSM.py -e 3 -l a -s 8 -o foo.png
EXPER = 1
SLIDING = 'b'
STEP = 1
figfilename = ""
ncfilename = ""
Mx = 601
try:
  opts, args = getopt(sys.argv[1:], "e:l:s:o:n:m:", \
                      ['exper=','sliding=','step=','out=','netcdf=','Mx='])
  for opt, arg in opts:
    if opt in ("-e", "--exper"):
      EXPER = float(arg)
    if opt in ("-l", "--sliding"):
      SLIDING = arg
    if opt in ("-s", "--step"):
      STEP = float(arg)
    if opt in ("-o", "--out"):
      figfilename = arg
    if opt in ("-n", "--netcdf"):
      ncfilename = arg
    if opt in ("-m", "--Mx"):
      Mx = int(arg)
except GetoptError:
  print "bad command line arguments ... exiting ..."
  sys.exit(-1)

SM_MODEL = 'B' # as defined in Schoof 2007, models A or B; ONLY USE B?

USE_SM_ERRORS = False

LL = 1800.0e3

#SM Specify your computational grid
user_grid = linspace(0.,LL,(Mx+1)/2);

#SM Parameters as defined in Schoof 2007, see also the MISMIP specifications
# for ALL MISMIP
n = 3.;              #SM Glen's law exponent
r = 0.9;            #SM ratio of ice to water density
rho_g = 8820.;        #SM 900 kg m^{-3} * 9.8 m s^{-2}
if USE_SM_ERRORS:
  a = 0.3 / (365. * 24. * 3600.); #### NON-STANDARD; SHOULD USE:
else:
  secpera = 31556926.  # (= 365.2422 days/a)
  a = 0.3 / secpera

if SLIDING == 'a':
  m = 1./3.;
  C = 7.624e6;      #SM m = 1/3 value
else:
  m = 1.;
  C = 7.2082e+010;  #SM m = 1 value

#SM Glen's law parameter in SI units
Alist = array(
        [[0.0, 4.6416e-24,  2.1544e-24,  1.0e-24,
               4.6416e-25,  2.1544e-25,  1.0e-25,
               4.6416e-26,  2.1544e-26,  1.0e-26,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
         [0.0, 3.0e-25, 2.5e-25, 2.0e-25,
               1.5e-25, 1.0e-25, 5.0e-26,
               2.5e-26, 5.0e-26, 1.0e-25,
               1.5e-25, 2.0e-25, 2.5e-25,
               3.0e-25,
          0.0, 0.0],
         [0.0, 1.6e-24, 1.4e-24, 1.2e-24,
               1.0e-24, 8.0e-25, 6.0e-25,
               4.0e-25, 2.0e-25, 4.0e-25,
               6.0e-25, 8.0e-25, 1.0e-24,
               1.2e-24, 1.4e-24, 1.6e-24]])   #  15th VALUE LABELED AS 16 IN Table 6 !?

if (EXPER == 1) | (EXPER == 2):
  A=Alist[0,STEP]
if (EXPER == 3) & (SLIDING == 'a'):
  A=Alist[1,STEP]
if (EXPER == 3) & (SLIDING == 'b'):
  A=Alist[2,STEP]

#SM models A and B, respectively, in Schoof 2007
if SM_MODEL == 'A':
  theta = 1.       #SM cold (theta = 1)
else:
  theta = 0.       #SM or warm (theta = 0) boundary layer specifications,

#print m,C,A,theta
#m = 1.; C = 7.2082e+010; A = 4.6416e-24; theta = 0.;


## from SMcold_bedheight.m:
def SMcold_bedheight(x):
  #SM DEPTH OF BED BELOW SEA LEVEL, i.e. gives b(x) in Schoof 2007.
  #SM NOTE the sign convention, SMcold_bedheight is positive if the bed is below
  #SM sea level, negative if above sea level.
  if USE_SM_ERRORS:
    #SM z =  -720 +778.5*(x/7e5);  ## SM CONTAINS TYPO!; should be "7.5e5"
    xx = x / 7.e5
  else:
    xx = x / 7.5e5
  if (EXPER == 3):
    return -(729. - 2184.8 * xx**2. + 1031.72 * xx**4. - 151.72 * xx**6.);
  else:
    return -720. + 778.5 * xx


## from SMcold_bedslope.m:
def SMcold_bedslope(x):
  #SM FIRST DERIVATIVE OF DEPTH OF BED BELOW SEA LEVEL; must agree with
  #SM SMcold_bedheight().
  if (EXPER == 3):
    xx = x / 7.5e5
    return -(-2184.8 * (2./7.5e5) * xx + 1031.72 * (4./7.5e5) * xx**3. \
             - 151.72 * (6./7.5e5) * xx**5.);
  else:
    return 778.5 / 7.5e5;


## from SMcold_function.m:
def SMcold_function(x):
  #SM Evaluates function whose zeros define x_g in `cold' steady marine sheet
  #SM problem
  h_f = r**(-1.) * SMcold_bedheight(x);
  b_x = SMcold_bedslope(x);
  s = a * x;
  return theta * a + C * s**(m+1.) / (rho_g * h_f**(m+2.)) - theta * s * b_x / h_f \
         - A * (rho_g * (1.-r)/4.)**n * h_f**(n+1.)


#SM Newton iteration parameters
#SM initial guess of grounding line position (metres)
if EXPER == 3:
  x = 800.0e3
else:
  x = 1270.0e3
delta_x = 10.;       #SM finite difference step size (metres) for gradient calculation
tolf = 1.e-4;        #SM tolerance for finding zeros
eps = finfo(float).eps    ### same as Matlab's "eps"
normf = tolf + eps;
toldelta = 1.e1;     #SM Newton step size tolerance
dx = toldelta + 1.;

#SM computes steady state grounding line position x based
#SM on intial guess in line 35, using equations (20) and
#SM (24) in Schoof 2007
print "Newton iteration for grounding line ...",
#print "  ",x,
while ((normf > tolf) | (abs(dx) > toldelta)):
    f = SMcold_function(x);
    normf = abs(f);
    grad = (SMcold_function(x+delta_x)-f)/delta_x;
    dx = -f/grad;
    x = x + dx;
    #print "  ",x,
print " x_g = %.3f km" % (x/1.e3)


#SM Calculate ice surface elevation S_Soln for grounded sheet

#SM grounded = user_grid(user_grid < x).';
grounded = transpose(user_grid[user_grid < x])

#SM specifies the grid points on which  ice surface elevations are to 
#SM be calculated; these are in reverse order, i.e. counting
#SM backwards from the grounding line at x to the ice
#SM divide at zero. FIRST POINT MUST BE GROUNDING LINE
#SM POSITION x COMPUTED ABOVE, EVEN IF THIS IS NOT
#SM INCLUDED IN YOUR COMPUTATIONAL GRID
#SM x_grid = [x; grounded(length(grounded):-1:1)];
x_grid = append([x],grounded[::-1])

#SM computes ice thickness at the grounding line as
#SM initial condition
h_f = SMcold_bedheight(x)/r;

def SMsurface(h,x):
  b_x = SMcold_bedslope(x)
  s = a * abs(x) # ???
  return b_x - (C / rho_g) * s**m / h**(m+1)

#SM computes ice THICKNESS H_soln at points with position X_soln
#options = odeset('AbsTol',1e-6,'RelTol',1e-6); #SM odeset('AbsTol',f/1e-3);
#[X_soln,H_soln] = ode45(@SMsurface,x_grid,[h_f],options);
#Hresult = integrate.odeint(SMsurface,[h_f],x_grid,atol=1.e-6,rtol=1.e-6)               
Hresult = integrate.odeint(SMsurface,[h_f],x_grid,atol=1.e-9,rtol=1.e-9)               
H_soln = array(Hresult[:,0])

#SM computes ice surface elevation by adding bed elevation
#S_soln = H_soln - SMcold_bedheight(X_soln);
S_soln = H_soln - SMcold_bedheight(x_grid);

#print "showing grounded part ..."
#plot(X_soln,-SMcold_bedheight(X_soln),'k')
#plot(X_soln,S_soln,'b')
plot(x_grid/1.e3,-SMcold_bedheight(x_grid),'k',linewidth=2)
hold(True)
plot(x_grid/1.e3,S_soln,'b',linewidth=2)
plot(x_grid/1.e3,zeros(size(x_grid)),'r:',linewidth=1) # show sea level

#SM H_soln is now rearranged and the point corresponding to the grounding
#SM line location is expunged to input these data into your computational grid:
#S_soln = S_soln(length(S_soln):-1:2);
#H_soln = H_soln(length(H_soln):-1:2);
S_soln = S_soln[:0:-1];
H_soln = H_soln[:0:-1];

#SM Calculate ice thickness for shelf from van der Veen (1986)
floating = transpose(user_grid[user_grid > x]);
q_0 = a * x;

#SM H_soln2 = h_f*(q_0 + a*(floating-x))./ ...
#SM    (q_0^(n+1) + h_f^(n+1)*((1-r)*rho_g/4)^n*A*((q_0 + a*(floating-x)).^(n+1) - q_0^(n+1))/a).^(1/(n+1)); 

numer = h_f * (q_0 + a * (floating - x))
base = q_0**(n+1.) + h_f**(n+1) * ((1.-r) * rho_g/4)**n * A \
                      * ((q_0 + a * (floating-x))**(n+1) - q_0**(n+1)) / a
H_soln2 = numer / (base**(1./(n+1.)));

print "ice thickness at x=0 km (divide): h_0 = %.3f m" % H_soln[0]
print "ice thickness at the grounding line: h_f = %.3f m" % h_f
if (x < 1800.0e3):
  print "ice thickness at x=1800km (calving front): h_c = %.3f m" % H_soln2[-1]
else:
  print "(grounding line past 1800km)"

user_thickness = append(H_soln,H_soln2)
#figure; plot(user_grid,user_thickness,'.-'); show()

x_gridfloat = append([x],floating)
thkfloat = append([h_f],H_soln2)
plot(x_gridfloat / 1.e3,(1.-r) * thkfloat,'b',linewidth=2);
plot(x_gridfloat / 1.e3,-r * thkfloat,'b',linewidth=2);
plot(x_gridfloat / 1.e3,-SMcold_bedheight(x_gridfloat),'k',linewidth=2) # show rest of bed
plot(x_gridfloat / 1.e3,zeros(size(x_gridfloat)),'r:',linewidth=1) # show sea level
hold(False)
xlabel("$y$  (km)",size=14)
ylabel("elevation  (m)",size=14)
axis([0.,1800.,-1000.,5000.])

# if .nc output *not* desired, either show figure or save it as PNG
if len(ncfilename) == 0:
  if len(figfilename) == 0:
    print "showing figure; close figure to end ..."
    show()
  elif figfilename == "0":
    print "(suppressing figure)"
  else:
    print "saving figure %s" % figfilename
    savefig(figfilename, dpi=300)
  sys.exit(0)


################## NetCDF write option #################
#### write a PISM-readable file with thickness only ####

from netCDF3 import Dataset as NC

My = 3
dy = (2. * LL) / (float(Mx-1))
Lx = (dy * float(My)) / 2.  # truely periodic in x direction

ncfile = NC(ncfilename, 'w')

# set global attributes
historysep = ' '
historystr = time.asctime() + ': ' + historysep.join(sys.argv) + '\n'
setattr(ncfile, 'history', historystr)
setattr(ncfile, 'source', "PISM: examples/mismip/solverMS.py")

# define the dimensions & variables
xdim = ncfile.createDimension('x', Mx)
ydim = ncfile.createDimension('y', My)
xvar = ncfile.createVariable('x', 'f8', dimensions=('x',))
yvar = ncfile.createVariable('y', 'f8', dimensions=('y',))
Hvar = ncfile.createVariable('thk', 'f4', dimensions=('y', 'x'))

# attributes of the variables
setattr(xvar, 'axis', 'X')
setattr(xvar, 'long_name', 'x-coordinate in Cartesian system')
setattr(xvar, 'standard_name', 'projection_x_coordinate')
setattr(xvar, 'units', 'm')

setattr(yvar, 'axis', 'Y')
setattr(yvar, 'long_name', 'y-coordinate in Cartesian system')
setattr(yvar, 'standard_name', 'projection_y_coordinate')
setattr(yvar, 'units', 'm')

setattr(Hvar, 'long_name', 'land ice thickness')
setattr(Hvar, 'standard_name', 'land_ice_thickness')
setattr(Hvar, 'units', 'm')

# write the dimension var data to the NetCDF file
for i in range(Mx):
  xvar[i]=-1800.0e3 + float(i) * dy
for i in range(My):
  yvar[i]=100.0* (-Lx + (float(i)+0.5) * dy)
  #xvar[i]=-Lx + float(i) * dy

# write the thickness data
thkrow = append(user_thickness[:0:-1],user_thickness)
thk = array([thkrow,thkrow,thkrow])  # uses My=3
Hvar[:] = thk

ncfile.close()
print "NetCDF file ",ncfilename," created with a single variable called 'thk'"

