#!/usr/bin/env python

from numpy import *
from matplotlib.pyplot import *
from sys import exit

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

nameroot = "lakes"

for dx in ("100","50","25"):
  basename = nameroot + dx + "km"
  filename = basename + ".nc"
  try:
    nc = NC(filename, 'r')
  except:
    print "can't read from file '%s' ... ENDING" % filename
    exit(-1)
  print "loading 'bwat' from %s ..." % filename

  xvar = nc.variables["x"]
  yvar = nc.variables["y"]
  x = asarray(squeeze(xvar[:]))
  y = asarray(squeeze(yvar[:]))

  bwatvar = nc.variables["bwat"]
  #print bwatvar
  bwat = asarray(squeeze(bwatvar[:])).transpose()
  print "  max = %9.3f, av = %8.3f" % (bwat.max(),bwat.sum()/(x.size*y.size))

  print "generating pcolor() image ..."
  pcolor(x/1000.0,y/1000.0,bwat,vmax=650.0)
  colorbar()
  #axis('equal')
  gca().set_aspect('equal')
  gca().autoscale(tight=True)
  xlabel('x  (km)')
  ylabel('y  (km)')

  nc.close()

  if len(dx)==3:
    pngfilename = "bwat_" + dx + "km" + ".png"
  else:
    pngfilename = "bwat_" + "0" + dx + "km" + ".png"
  print "saving figure in %s ..." % pngfilename
  savefig(pngfilename, dpi=300, bbox_inches='tight')
  #show()
  close()

