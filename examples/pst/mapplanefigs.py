#!/usr/bin/env python

from pylab import *
import pycdf
import os
import sys

def savepng(fignum,prefix):
  filename = prefix + ".png"
  print "saving figure(%d) as %s ..." % (fignum,filename)
  savefig(filename, dpi=300, facecolor='w', edgecolor='w')

for j in [1,2,3,4]:
  filename = "P%dfinest.nc" % j
  print "opening %s to make map plane csurf figure" % filename
  nc = pycdf.CDF(filename)
  x = nc.var("x").get()
  y = nc.var("y").get()
  csurf = squeeze(nc.var("csurf").get())
  print "  csurf has max = %.2f and min = %.2f (m/a)" % (csurf.max(),csurf.min())
  nc.close()

  # show speed as color
  figure(j, figsize=(9,8))
  #figure(j, figsize=(9,8));clf();hold(True)
  if (j == 1):
    pcolor(x,y,log10(maximum(csurf,1.0)),shading='flat',
           vmin=log10(5.0),vmax=log10(300.0))
    commonCmap = cm.get_cmap()
  else:
    pcolor(x,y,log10(maximum(csurf,1.0)),shading='flat',
           vmin=log10(5.0),vmax=log10(300.0),
           cmap = commonCmap)
  cb = colorbar(ticks = log10([10,25,50,100,250]),
                #ticks=[0,log10(5),1,log10(50),2,log10(500.0)],
                aspect=15,shrink=0.85, extend='both',
                format=FormatStrFormatter('$10^{%.3f}$'))
  axis('off')
  axis('equal')

  prefix = 'P%d_5km_csurf' % j
  savepng(j,prefix)
  print "autocropping %s.png" % prefix
  # uses one of the ImageMagick tools (http://www.imagemagick.org/)
  try:
    status = os.system("mogrify -verbose -trim +repage %s.png" % prefix)
  except KeyboardInterrupt:  sys.exit(2)
  if status:  sys.exit(status)


