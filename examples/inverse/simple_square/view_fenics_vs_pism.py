import matplotlib.pyplot as pp
import numpy as np
from Scientific.IO.NetCDF import NetCDFFile as Dataset

ncfile_pism = Dataset('simple_square.nc','r') 
ncfile_fenics = Dataset('/home/hmarijke/UAF-inverse/maxwell_inverse/python/experiments/simple_square/results/square_fine_comp_pism.nc','r')

# read out tauc for pism run and gamma for fenics run
tauc = ncfile_pism.variables['tauc'][:]
gamma = ncfile_fenics.variables['gamma_IGN'][:]

# # manipulate tauc so that it's the same shape and orientation as gamma
tauc = tauc[0,:,:]
tauc = tauc.transpose()
tauc = tauc[:,::-1]

# get the same units for both measurements 
# gamma is in bar a/m
# tauc is in Pa s/m
tauc = tauc * 1.e-5 / 31556925.9747 # 31556925.9747 sec per year

# plot the difference between the two
pp.figure()
pp.imshow((gamma-tauc)/gamma.max()*100, origin='lower')
pp.colorbar()
pp.title('perc. diff. of \n gamma and tauc')


maxG = 0.0055
minG = 0.0001
# plot tauc
pp.figure()
pp.imshow(tauc, origin='lower') #, vmin=minG,vmax=maxG)
pp.colorbar()
pp.title('tauc \n bar a/m')


# plot gamma
pp.figure()
pp.imshow(gamma, origin='lower', vmin=minG,vmax=maxG)
pp.colorbar()
pp.title('gamma \n bar a/m')


###############################
## compare the two u_true
# read out u_true for pism and fenics run
ux_true_pism = ncfile_pism.variables['u_true'][:]
uy_true_pism = ncfile_pism.variables['v_true'][:]
ux_true_pism = ux_true_pism[0,:,:]
ux_true_pism = ux_true_pism.transpose()
uy_true_pism = uy_true_pism[0,:,:]
uy_true_pism = uy_true_pism.transpose()
u_true_pism = np.zeros((121,61))
for i in range(121):
    for j in range(61):
        u_true_pism[i,j] = np.sqrt(ux_true_pism[i,j]**2 + uy_true_pism[i,j]**2)
u_true_pism = u_true_pism[:,::-1]
#u_true_pism = u_true_pism * 31556925.9747 # PISM vel is in m/s -> change to m/a

ux_true_fenics = ncfile_fenics.variables['ux_true'][:]
uy_true_fenics = ncfile_fenics.variables['uy_true'][:]
u_true_fenics = np.zeros((121,61))
for i in range(121):
    for j in range(61):
        u_true_fenics[i,j] = np.sqrt(ux_true_fenics[i,j]**2 + uy_true_fenics[i,j]**2)

# calculate the difference and plot
diff = u_true_fenics - u_true_pism
percent_diff = diff / u_true_fenics.max() *100
pp.figure()
pp.imshow(percent_diff, origin='lower')
pp.colorbar()
pp.title('perc. diff. of \n u_true_fenics and u_true_pism')


pp.figure()
pp.imshow(u_true_fenics, origin='lower')
pp.colorbar()
pp.title('u_true_fenics \n m/a')


pp.figure()
pp.imshow(u_true_pism, origin='lower')
pp.colorbar()
pp.title('u_true_pism \n m/a')



###############################
## compare the two tauc_true
# read out gamma_true and tauc_true for pism and fenics run
tauc_true_pism = ncfile_pism.variables['tauc_true'][:]
tauc_true_pism = tauc_true_pism[0,:,:]
tauc_true_pism = tauc_true_pism.transpose()
tauc_true_pism = tauc_true_pism[:,::-1]
tauc_true_pism = tauc_true_pism * 1.e-5 / 31556925.9747

gamma_true_fenics = ncfile_fenics.variables['gamma_true'][:]


# calculate the difference and plot
diff = gamma_true_fenics - tauc_true_pism
pp.figure()
pp.imshow(diff, origin='lower')
pp.colorbar()
pp.title('gamma_true_fenics - tauc_true_pism \n bar a/m')


# plot tauc_true
pp.figure()
pp.imshow(tauc_true_pism,origin='lower')
pp.colorbar()
pp.title('tauc_true pism \n bar a/m')

# plot gamma_true
pp.figure()
pp.imshow(gamma_true_fenics,origin='lower')
pp.colorbar()
pp.title('gamma_true fenics \n bar a/m')
pp.show()
