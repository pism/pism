#!/usr/bin/env python

from pylab import *
import pycdf
import os
import sys

prefixlist = []

def savepng(fignum,prefix):
  prefixlist.append(prefix)
  filename = prefix + ".png"
  print "saving figure(%d) as %s ..." % (fignum,filename)
  savefig(filename, dpi=300, facecolor='w', edgecolor='w')

## csurf for each experiments
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

## bwat for P4
filename = "P4finest.nc"
print "opening %s to make map plane bwat figure" % filename
nc = pycdf.CDF(filename)
x = nc.var("x").get()
y = nc.var("y").get()
bwat = squeeze(nc.var("bwat").get())
print "  bwat has max = %.2f and min = %.2f (m/a)" % (bwat.max(),bwat.min())
nc.close()
figure(5, figsize=(9,8))
pcolor(x,y,bwat,shading='flat',vmin=0.0,vmax=2.0)
colorbar(ticks = [0.0,0.5,1.0,1.5,2.0],
         aspect=15,shrink=0.85, extend='both')
axis('off')
axis('equal')
prefix = 'P4_5km_bwat'
savepng(5,prefix)

## csurf for P4 at 5 years
filename = "P4finest_5yr.nc"
print "opening %s to make map plane csurf figure" % filename
nc = pycdf.CDF(filename)
x = nc.var("x").get()
y = nc.var("y").get()
csurf = squeeze(nc.var("csurf").get())
print "  csurf has max = %.2f and min = %.2f (m/a)" % (csurf.max(),csurf.min())
nc.close()
figure(6, figsize=(9,8))
pcolor(x,y,log10(maximum(csurf,1.0)),shading='flat',
           vmin=log10(5.0),vmax=log10(3600.0))
cb = colorbar(ticks = log10([10,30,100,300,1000,3000]),
              aspect=15,shrink=0.85, extend='both',
              format=FormatStrFormatter('$10^{%.3f}$'))
axis('off')
axis('equal')
savepng(6,'P4_5km_5yr_csurf')

## csurf detail for P1, stream 0 (on 5km grid)
filename = "P1finest.nc"
print "opening %s to make map plane csurf detail figure" % filename
nc = pycdf.CDF(filename)
x = nc.var("x").get()
y = nc.var("y").get()
csurf = squeeze(nc.var("csurf").get())
xpatch = x[200:281] # 301 x 301
ypatch = y[130:171]
csurfpatch = csurf[130:171,200:281]
print "  csurfpatch has max = %.2f and min = %.2f (m/a)" % (csurfpatch.max(),csurfpatch.min())
nc.close()
figure(7, figsize=(9,8))
pcolor(xpatch,ypatch,log10(maximum(csurfpatch,1.0)),shading='flat',
           vmin=log10(5.0),vmax=log10(300.0))
# use colorbar from whole csurf above
#cb = colorbar(ticks = log10([10,25,50,100,250]),
#              aspect=15,shrink=0.85, extend='both',
#              format=FormatStrFormatter('$10^{%.3f}$'))
axis('off')
axis('equal')
savepng(7,'P1detail_5km')

## csurf detail for P1, stream 0 (on 7.5km grid)
filename = "P1fine.nc"
print "opening %s to make map plane csurf detail figure" % filename
nc = pycdf.CDF(filename)
x = nc.var("x").get()
y = nc.var("y").get()
csurf = squeeze(nc.var("csurf").get())
xpatch = x[134:188] # 201 x 201
ypatch = y[87:114]
csurfpatch = csurf[87:114,134:188]
print "  csurfpatch has max = %.2f and min = %.2f (m/a)" % (csurfpatch.max(),csurfpatch.min())
nc.close()
figure(8, figsize=(9,8))
pcolor(xpatch,ypatch,log10(maximum(csurfpatch,1.0)),shading='flat',
           vmin=log10(5.0),vmax=log10(300.0))
axis('off')
axis('equal')
savepng(8,'P1detail_7.5km')

## csurf detail for P1, stream 0 (on 10km grid)
filename = "P1.nc"
print "opening %s to make map plane csurf detail figure" % filename
nc = pycdf.CDF(filename)
x = nc.var("x").get()
y = nc.var("y").get()
csurf = squeeze(nc.var("csurf").get())
xpatch = x[100:141] # 151 x 151
ypatch = y[65:86]
csurfpatch = csurf[65:86,100:141]
print "  csurfpatch has max = %.2f and min = %.2f (m/a)" % (csurfpatch.max(),csurfpatch.min())
nc.close()
figure(9, figsize=(9,8))
pcolor(xpatch,ypatch,log10(maximum(csurfpatch,1.0)),shading='flat',
           vmin=log10(5.0),vmax=log10(300.0))
axis('off')
axis('equal')
savepng(9,'P1detail_10km')

## csurf detail for P1, stream 0 (on 15km grid)
filename = "P1coarse.nc"
print "opening %s to make map plane csurf detail figure" % filename
nc = pycdf.CDF(filename)
x = nc.var("x").get()
y = nc.var("y").get()
csurf = squeeze(nc.var("csurf").get())
xpatch = x[66:95] # 101 x 101
ypatch = y[43:58]
csurfpatch = csurf[43:58,66:95]
print "  csurfpatch has max = %.2f and min = %.2f (m/a)" % (csurfpatch.max(),csurfpatch.min())
nc.close()
figure(10, figsize=(9,8))
pcolor(xpatch,ypatch,log10(maximum(csurfpatch,1.0)),shading='flat',
           vmin=log10(5.0),vmax=log10(300.0))
axis('off')
axis('equal')
savepng(10,'P1detail_15km')


## surface elevation profile for P1_100k versus P0A [NOT ACTUALLY MAP-PLANE]
filename = "P1_100k.nc"
print "opening %s to make surface elevation profile figure" % filename
nc = pycdf.CDF(filename)
x = nc.var("x").get()
y = nc.var("y").get()
usurf = squeeze(nc.var("usurf").get())
print "  usurf has max = %.2f and min = %.2f (m) and shape =" % (usurf.max(),usurf.min()),
print shape(usurf)
nc.close()
filename = "P0A.nc"
print "opening %s to make surface elevation profile figure" % filename
nc = pycdf.CDF(filename)
x0 = nc.var("x").get()
y0 = nc.var("y").get()
usurf0 = squeeze(nc.var("usurf").get())
print "  usurf0 has max = %.2f and min = %.2f (m) and shape =" % (usurf0.max(),usurf0.min()),
print shape(usurf0)
nc.close()
figure(99, figsize=(9,6))
mylabel = "P1"
plot(x[75:151]/1000.0,usurf[75,75:151]-2000.0,'k.-',
     linewidth=1.5,label=mylabel,markersize=10.)
hold(True)
mylabel = "initial state"
plot(x0[75:151]/1000.0,usurf0[75,75:151]-2000.0,'k:',
     linewidth=2.5,label=mylabel)
hold(False)
legend(loc='upper right')
axis([0, 750, 0, 4000]) 
xticks(arange(0,800,100),fontsize=14)
yticks(arange(0,5000,1000),fontsize=14)
xlabel("x  (km)",fontsize=16)
ylabel("surface elevation  (m)",fontsize=16)
savepng(99,'profile_P1_P0A')


print ""
print "AUTOCROPPING with mogrify"
# uses one of the ImageMagick tools (http://www.imagemagick.org/)
#for filepre in ["P1_5km_csurf","P2_5km_csurf","P3_5km_csurf","P4_5km_csurf",\
#                "P4_5km_bwat",]:
for filepre in prefixlist:
  try:
    status = os.system("mogrify -verbose -trim +repage %s.png" % filepre)
  except KeyboardInterrupt:  sys.exit(2)
  if status:  sys.exit(status)


