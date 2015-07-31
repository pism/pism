#!/usr/bin/env python

#!/usr/bin/env python

import netCDF4 as nc
import sys
import matplotlib.pylab as mpl
import numpy as np
import seaborn as sns

def main ():
    tracers=nc.Dataset("tracers_traced.nc", "r")
    tx = tracers.variables["x"][:,0]
    ty = tracers.variables["y"][:,0]
    tz = tracers.variables["z"][:,0]
    tt = tracers.variables["time"][:]
    sf = nc.Dataset("stress_file.nc", "r")
    sx = sf.variables["x"][:]
    sy = sf.variables["y"][:]
    index = 0
    for (n, x) in enumerate(sx):
        if x == 3e5:
            index = n
            break
    if (index == 0 ):
        print ("Unuseable input file. x does not have a value 3e5\n")
        sys.exit(666)
    # slope at index should be 0.02 in y dir.
    if (sy[len(sy)/2] !=0 ):
        print ('Unuseable input file. Must have y=0 at [len/2]')
        sys.exit(666)
    if (np.abs(sf.variables["h_y_j"][0,index,len(sy)/2] / 0.02 -1 ) > 0.01):
        print ("trouble with slope! should be 0.02 at (%d, %d)  "%(index, len(sy)/2))
        print ("slope is %f"%(sf.variables["h_y_j"][0,index,len(sy)/2]))
        sys.exit(666)
        
        
    phi_obs = np.arctan2(ty,tx)
    
    A = 3.1689e-24
    secpera = 365*86400
    phi_comp = -(.5* (910* 9.81 *0.02/3e5)**3 * A * 3e5**2 * (1000**4 - 0.2 *((0.109891/365.2425/86400. * tt )**4))*tt)%(2*np.pi)
    phi_comp = phi_comp - 2 * np.pi * (phi_comp > np.pi)
    mpl.figure()
    mpl.plot(tt/secpera, phi_obs-phi_comp , '.', label="diff")
    mpl.xlabel("time in years")
    mpl.ylabel("$\\varphi_{modeled} - \\varphi_{expect}$")
    mpl.title("Error in angle along circle")
    mpl.savefig("E_phi.png")
    mpl.figure()
    mpl.plot(tt/secpera, phi_obs, '.', label = "model")
    mpl.plot(tt/secpera, phi_comp, '.', label = "analytic")
    mpl.xlabel("time in years")
    mpl.ylabel("$\\varphi$")
    mpl.title("Angle along circle")
    mpl.legend()
    mpl.savefig("phi.png")
    mpl.figure()
    mpl.plot(tt/secpera, np.sqrt(tx**2+ty**2)-3e5)
    mpl.xlabel("time in years")
    mpl.ylabel("Error in distance from center (R=300000 m)")
    mpl.title("Circle radius error")
    mpl.savefig("E_r.png")
    mpl.figure()
    mpl.plot(tt/secpera, tz - (1000 - 0.1/.91/365.2425/86400. * tt ))
    mpl.xlabel("time in years")
    mpl.ylabel("Error in vertical position (m)")
    mpl.title("Vertical error")
    mpl.savefig("E_z.png")
    mpl.show()
            #    if (sf.variables[""])
            
            # stress in x dir is const in x, only depends on y
            # stress in y dir is const in y, only depends on x
            # therefore: get y stress at j point and x stress at i point (the others are shifted)
            
if __name__ == "__main__":
    
    main()
    
