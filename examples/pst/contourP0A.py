#!/usr/bin/env python

# draw contour map of P0A thickness and image of P0A basal water thickness
# see Bueler and Brown (2008 preprint), "The shallow shelf approximation
#    as a sliding law ...", http://arxiv.org/abs/0810.3449
# optionally produces image files  thkP0Acontour.png, bwatP0Aimage.png

from pylab import *
from netCDF3 import Dataset as NC

filename = "P0A.nc"
print "opening %s to make contour map of thickness and image of bwat" % filename
nc = NC(filename, 'r')
x = nc.variables["x"][:]
y = nc.variables["y"][:]
thk = squeeze(nc.variables["thk"][:])
print "  thk has max = %.2f and min = %.2f (m/a)" % (thk.max(),thk.min())
bwat = squeeze(nc.variables["bwat"][:])
print "  bwat has max = %.2f and min = %.2f (m/a)" % (bwat.max(),bwat.min())
nc.close()

figure(1,figsize=(6,5))
thkContours = [1,500,1000,1500,2000,2500,3000,3500,3500]
thkCS = contour(x/1.0e3,y/1.0e3,thk,thkContours,colors='k',linewidths=1.5)
clabel(thkCS,[1000,2000,3000],fmt='%d')
axis('equal')
xmax=700
axis([-xmax,xmax,-xmax,xmax])
xlabel('x  (km)',fontsize=16)
ylabel('y  (km)',fontsize=16)
xticks([-500,0,500],fontsize=14)
yticks([-500,0,500],fontsize=14)

figure(2)
imshow(bwat,aspect='equal',interpolation='nearest',vmin=0,vmax=2.2)
axis('off')

print "close figures to finish"
show()

## optional save to PNG
#print "saving figure(1) as thkP0Acontour.png ..."
#savefig("thkP0Acontour.png", dpi=300, facecolor='w', edgecolor='w')
#print "saving figure(2) as bwatP0Aimage.png ..."
#savefig("bwatP0Aimage.png", dpi=300)
#print "remember to 'mogrify -trim +repage *.png'"

