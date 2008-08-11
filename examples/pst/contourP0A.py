#!/usr/bin/env python

# produces figures
#   thkP0Acontour.png
#   bwatP0Aimage.png

from pycdf import *
from pylab import *

nc = CDF("P0A.nc")
x = nc.var("x").get()
y = nc.var("y").get()
thk = squeeze(nc.var("thk").get())
bwat = squeeze(nc.var("bwat").get())

figure(1,figsize=(6,5))
#subplot(122)
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
savefig("thkP0Acontour.png", dpi=300, facecolor='w', edgecolor='w')
print "remember to 'mogrify -trim +repage thkP0Acontour.png'"

figure(2)
imshow(bwat,aspect='equal',interpolation='nearest',vmin=0,vmax=2.2)
axis('off')
savefig("bwatP0Aimage.png", dpi=300)
print "remember to 'mogrify -trim +repage bwatP0Aimage.png"

#show()

