#!/usr/bin/env python

from scipy import *
from pylab import *
import sys

# solverMISMIP.py is matlab --> pylab translation of "SMsolver2.m" and other
# Matlab codes contained in
#    http://homepages.ulb.ac.be/~fpattyn/mismip/MISMIP_distribution.tar
# see major comments at the start of SMsolver2.m; comments here 
# starting "#SM" are from the Matlab codes


#SM Specify your computational grid; user grid must be a COLUMN vector in metres
user_grid = linspace(0,1.8e6,150);
#user_grid = linspace(0,1.8e6,1500);


#SM Parameters as defined in Schoof 2007, see also the MISMIP specifications
n = 3.;              #SM Glen's law exponent
m = 1.;#SM1/3;         #SM sliding exponent
r = 0.9;            #SM ratio of ice to water density
rho_g = 8820.;        #SM 900 kg m^{-3} * 9.8 m s^{-2}
A = 4.6416e-24;               #SM Glen's law parameter in SI units
                   #SM NB Paterson (1994) does not use SI units!
#SMC = 7.624e6;        #SM m = 1/3 value
C = 7.2082e+010;     #SM m = 1 value 
                    #SM for both, tau_b = 80kPa gives u ~ 35 m a^(-1)
#a = 0.3 / (365. * 24. * 3600.);  
#### SHOULD USE
secpera = 31556926. # (= 365.2422 days/a)
a = 0.3 / secpera
theta = 0.;          #SM cold (theta = 1) or warm (theta = 0) boundary layer specifications,
                    #SM models A and B respectively in Schoof 2007
x = 1.270e6;        #SM initial guess of grounding line position (metres)
eps = finfo(float).eps    ### same as Matlab's "eps"

#SM Newton iteration parameters
delta_x = 10.;       #SM finite difference step size (metres) for gradient calculation
tolf = 1.e-4;        #SM tolerance for finding zeros
normf = tolf + eps;
toldelta = 1.e1;     #SM Newton step size tolerance
dx = toldelta + 1.;


## from SMcold_bedheight.m:
def SMcold_bedheight(x):
  #SM DEPTH OF BED BELOW SEA LEVEL, i.e. gives b(x) in Schoof 2007.
  #SM NOTE the sign convention, SMcold_bedheight is positive if the bed is below
  #SM sea level, negative if above sea level.
  #SM z =  -720 +778.5*(x/7e5);  ## SM CONTAINS TYPO!; should be "7.5e5"; TYPO! TYPO! TYPO! TYPO!
  #xx = x / 7.e5
  #### SHOULD USE
  xx = x / 7.5e5
  return -720. + 778.5 * xx
  #SM -(729 - 2184.8*(x/7.5e5).^2 + 1031.72*(x/7.5e5).^4 - 151.72*(x/7.5e5).^6);
  #return -(729. - 2184.8 * xx**2. + 1031.72 * xx**4. - 151.72 * xx**6.);


## from SMcold_bedslope.m:
def SMcold_bedslope(x):
  #SM FIRST DERIVATIVE OF DEPTH OF BED BELOW SEA LEVEL; must agree with
  #SM SMcold_bedheight.
  return 778.5 / 7.5e5;
  #SM -(-2184.8*2*x/7.5e5^2 + 1031.72*4*x.^3/7.5e5^4 - 151.72*6*x.^5/7.5e5^6);
  #xx = x / 7.5e5
  #return -(-2184.8 * (2./7.5e5) * xx + 1031.72 * (4./7.5e5) * xx**3. \
  #         - 151.72 * (6./7.5e5) * xx**5.);


## from SMcold_function.m:
def SMcold_function(x):
  #SM Evaluates function whose zeros define x_g in `cold' steady marine sheet
  #SM problem
  #global m n A C rho_g r a theta
  h_f = r**(-1.) * SMcold_bedheight(x);
  b_x = SMcold_bedslope(x);
  s = a * x;
  #z = theta*a + C*s.^(m+1)./(rho_g*h_f.^(m+2)) - theta*s.*b_x./h_f - A*(rho_g*(1-r)/4)^n*h_f.^(n+1);
  return theta * a + C * s**(m+1.) / (rho_g * h_f**(m+2.)) - theta * s * b_x / h_f \
         - A * (rho_g * (1.-r)/4.)**n * h_f**(n+1.)


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
print "ice thickness at the grounding line: h_f = %.3f m" % h_f

def SMsurface(h,x):
  b_x = SMcold_bedslope(x);
  s = a * x;
  return b_x - (C / rho_g) * s**m / h**(m+1);

#SM computes ice THICKNESS H_soln at points with
#SM position X_soln
#SM options = odeset('AbsTol',1e-6,'RelTol',1e-6); #SM odeset('AbsTol',f/1e-3);
#[X_soln,H_soln] = ode45(@SMsurface,x_grid,[h_f],options);
Hresult = integrate.odeint(SMsurface,[h_f],x_grid,atol=1.e-6,rtol=1.e-6)               
H_soln = array(Hresult[:,0])

#SM computes ice surface elevation by adding bed elevation
#S_soln = H_soln - SMcold_bedheight(X_soln);
S_soln = H_soln - SMcold_bedheight(x_grid);

#print "showing grounded part ..."
#plot(X_soln,-SMcold_bedheight(X_soln),'k')
#plot(X_soln,S_soln,'b')
plot(x_grid/1.e3,-SMcold_bedheight(x_grid),'k')
hold(True)
plot(x_grid/1.e3,S_soln,'b')
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

#print "showing floating part ..."
#plot([x; floating],(1-r)*[h_f; H_soln2],'b');
#plot([x; floating],-r*[h_f; H_soln2],'b');
x_gridfloat = append([x],floating) / 1.e3
thkfloat = append([h_f],H_soln2)
plot(x_gridfloat,(1.-r) * thkfloat,'b');
plot(x_gridfloat,-r * thkfloat,'b');
plot(x_gridfloat,zeros(size(x_gridfloat)),'r:',linewidth=1) # show sea level
hold(False)
xlabel("x  (m)")
ylabel("elevation  (m)")
#show()
outputname = "SMout.png"
print "saving figure %s" % outputname
savefig(outputname, dpi=300)

user_thickness = append(H_soln,H_soln2)
#figure; plot(user_grid,user_thickness,'.-'); show()


