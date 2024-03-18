#!/usr/bin/env python3

# Helper scripts for postprocessing PISM result to upload for CalvinMIP, 
# by albrecht@pik-potsdam.de

import numpy as np
import netCDF4 as nc


###############################################################################################

# define piecewise linear transects between i-j-locations 
def get_troughs(of,trans,dp):

  try:
    compfile = nc.Dataset(of, 'r')
  except Exception:
    print("Specify NetCDF input file.")
    print(of)
    exit(-1)
  xobs = np.squeeze(compfile.variables["x"][:])
  yobs = np.squeeze(compfile.variables["y"][:])
  compfile.close()

  trans_points=[]
  
  for l,points in enumerate(trans):

    pkm=[]
    for pp in points:
      pkm.append([xobs[pp[0]],yobs[pp[1]]])

    #dp=np.float(dkm)/np.float(resolution)
    dpr=dp

    po=[]
    po.append(points[0])
    pox=po[0][0]
    poy=po[0][1]

    count=0
    d=0

    while count<len(points)-1:

      xdir=np.sign(points[count+1][0]-points[count][0])
      ydir=np.sign(points[count+1][1]-points[count][1])
      my=np.float(points[count][1]-points[count+1][1])
      mx=np.float(points[count][0]-points[count+1][0])

      if mx==0:
        dpx=0.0
        dpy=ydir*dpr
        mnew=np.nan
      elif my==0:
        dpy=0.0
        dpx=xdir*dpr
        mnew=0.0
      else:  
        mnew=my/mx
        dpy=ydir*np.sqrt((mnew*dpr)**2/(1.0+mnew**2))
        dpx=dpy/mnew

      dx=dpx*(xobs[points[count][0]+xdir]-xobs[points[count][0]])
      dy=dpy*(yobs[points[count][1]+ydir]-yobs[points[count][1]])

      px=dpx+pox
      py=dpy+poy

      if (xdir*(px-points[count+1][0])<=0.0 and ydir*(py-points[count+1][1])<=0.0):
        po.append([px,py])
        dpr=dp
        pox=po[d+1][0]
        poy=po[d+1][1]
      else: #overshoot, next point, take rest distance
        dpr=np.sqrt((px-points[count+1][0])**2+(py-points[count+1][1])**2)
        pox=points[count+1][0]
        poy=points[count+1][1]
        count+=1
        d-=1
      d+=1

    trans_points.append(po)
  return trans_points


def interpolate_along_transect(var,i,j,di,dj): 
    return (1.0-dj)*((1.0-di)*var[i,j] + di*var[i,j+1]) + dj*((1.0-di)*var[i+1,j] + di*var[i+1,j+1])

def nearest_along_transect(var,i,j,di,dj):
    return var[i+np.int(np.round(di)),j+np.int(np.round(dj))]


