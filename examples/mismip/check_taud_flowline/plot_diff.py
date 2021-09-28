
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import netCDF4 as nc

# ICE FRONT 

#data_path = "/home/albrecht/pism21/pism-dev/examples/mismip/check_taud_flowline2/"
data_path = ""

experiments = [data_path+'result_1a_A1_800_300_1750e3_601_onesided_',
               data_path+'result_1a_A1_800_300_1750e3_601_centered_']
labels = ['one sided', 'centered']


# FIXME resolution dependent!!

refdat = data_path+'ICESHELF_1a_A1_800_300_1750e3_601_default.nc'
ncr    = nc.Dataset(refdat, 'r')
velref = np.ma.array(ncr.variables["u_ssa_bc"][:])
ncr.close()

ncr    = nc.Dataset(experiments[0]+'default.nc', 'r')
mask = np.ma.array(ncr.variables["mask"][0,:])
x = np.ma.array(ncr.variables["x"][:])
#thk_ref = np.ma.array(ncr.variables['thk'][0,:])
#print(np.shape(mask))
ncr.close()

divide_idx = int(np.floor(mask.shape[1]/2))
divide_idx = 2

glidxs = np.where(mask[1,divide_idx:]==3)
glidx = glidxs[0][0]+divide_idx -1 # Last grounded cell

ifidxs = np.where(mask[1,divide_idx:]==4)
ifidx = ifidxs[0][0]+divide_idx -1


#MISMIP
secpera=3.15569259747e7


fig, axes = plt.subplots(4, 2, sharex=False, sharey=False, figsize=(10, 8))
plt.subplots_adjust( wspace=0.25, hspace=0.1)

for i,exp in enumerate([ 'default.nc', 'p1.nc', 'p2.nc', 'p3.nc']):
    #print(i, 2*i, 2*i+1,exp)
    # one figure with changes in ssa vels and one with taud
    ax2 = axes.flatten()[2*i]
    ax1 = axes.flatten()[2*i+1]
    #experiments = glob.glob(os.path.join(data_path, ensemble_id+"_gradienthaseloff/diagnostic_*_"+exp))
    varname = "u_ssa"
    
    ax1.axvline(x=x[glidx]/1000., color='grey', linestyle='--', label=None)
    ax1.axhline(y=0, color='grey', linestyle='--', label=None)
    for n,run in enumerate(experiments):
        refvariable=velref

        ncr    = nc.Dataset(run+exp, 'r')
        variable = np.ma.array(ncr.variables[varname][0,:])*secpera
        thk = np.ma.array(ncr.variables['thk'][0,:])
        ncr.close()

        variable[mask==4] = np.nan
        ax1.plot(x[divide_idx:]/1000, variable[1,divide_idx:]-refvariable[1,divide_idx:],
                 linewidth=2*(2-n),
                 label=labels[n])

    ax1.set_ylim([-2000,500]) # km
    ax1.set_xlim([0,1760]) # km

    
    
    varname = "taud_x"
    #varname = "taud_mag"
    ax2.axvline(x=(x[glidx]), color='grey', linestyle='--', label=None)
    ax2.axhline(y=0, color='grey', linestyle='--', label=None)      
    for n,run in enumerate(experiments):
        ncr    = nc.Dataset(run+"default.nc", 'r')
        refvariable = np.ma.array(ncr.variables[varname][0,:])
        ncr.close()

        ncr    = nc.Dataset(run+exp, 'r')
        variable = np.ma.array(ncr.variables[varname][0,:])
        ncr.close()


        variable[mask==4] = np.nan
        ax2.plot(x[divide_idx:]/1000, variable[1,divide_idx:]-refvariable[1,divide_idx:],
                    linewidth=2*(2-n),
                    label=labels[n])
        ax2.scatter(x[divide_idx:]/1000, variable[1,divide_idx:]-refvariable[1,divide_idx:], label=None)
    
    if i==1:    
      ax2.set_ylabel('Change in driving stress (Pa)')   
      ax1.set_ylabel('Change in velocity (m/a)')
    ax2.set_ylim([-1200,600]) # km
    ax2.set_xlim([1700,1755]) # km      
    ax2.text(1702,400,exp.split('_')[-1].rstrip('.nc'))
    if i<3: 
      ax1.set_xticklabels([])
      ax2.set_xticklabels([]) 
 
ax = axes.flatten()[0]
ax.legend(loc='lower right')
axes.flatten()[6].set_xlabel('Distance from ice divide (km)')
axes.flatten()[7].set_xlabel('Distance from ice divide (km)')

#plt.savefig('/home/albrecht/www/pism_taudfix/taudfix_margin.png', bbox_inches='tight')
plt.show()
