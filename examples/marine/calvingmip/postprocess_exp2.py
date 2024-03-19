#!/usr/bin/env python3

# Creates output file from PISM EXP2 result to upload for CalvinMIP, 
# as instructed from https://github.com/JRowanJordan/CalvingMIP/wiki/Experiment-2

import numpy as np
import netCDF4 as nc
import datetime
import postprocess_helper as ph

secperyear=365*24*3600

resolution=5.0
dkm=5.0 #km steps

pismpath     = "/p/tmp/albrecht/pism23/calvmip/circular/exp2-05km-dir/"
pism_outfile = pismpath + "results/extra_exp2.nc"
pism_tsfile  = pismpath + "results/ts_exp2.nc"

pism_infile = pismpath + "input/circular_input_5km.nc"

exp_outfile = "CalvingMIP_EXP2_PISM_PIK.nc"

#################################################################################################

print('Load PISM data...')
datnc = nc.Dataset(pism_outfile,"r")
exp_x = datnc.variables["x"][:]
exp_y = datnc.variables["y"][:]
exp_t = datnc.variables["time"][:]/secperyear-120000.0
try:
    exp_xvm = datnc.variables["xvelmean"][:]*secperyear
    exp_yvm = datnc.variables["yvelmean"][:]*secperyear
    exp_thk = datnc.variables["lithk"][:]
except:
    exp_xvm = datnc.variables["u_ssa"][:]*secperyear
    exp_yvm = datnc.variables["v_ssa"][:]*secperyear
    exp_thk = datnc.variables["thk"][:]
exp_mask = datnc.variables["mask"][:]-1
exp_crate = datnc.variables["calvingmip_calving_rate"][:]*secperyear
exp_topg = datnc.variables["topg"][:]
datnc.close()

#print(exp_t)
Mt,Mx,My = np.shape(exp_mask)

exp_thk[exp_thk == 0.0] = np.nan
exp_xvm[exp_xvm == 0.0] = np.nan
exp_yvm[exp_yvm == 0.0] = np.nan

datnc = nc.Dataset(pism_tsfile,"r")
exp_ts_t      = [datnc.variables["time"][:]/secperyear-120000.0]
exp_ts_afl    = [datnc.variables["iareafl"][:]]
exp_ts_agr    = [datnc.variables["iareagr"][:]]
exp_ts_lim    = [datnc.variables["lim"][:]]
exp_ts_limnsw = [datnc.variables["limnsw"][:]]
exp_ts_tendcf = [datnc.variables["tendlicalvf"][:]*secperyear]
exp_ts_tendgf = [datnc.variables["tendligroundf"][:]*secperyear]
datnc.close()
#print(exp_ts_t)

###############################################################################################

# define profiles
Mp=int((Mx-1)/2.0)
print(Mx,My,Mp)

points_C = [[Mp,Mp],[Mx-2,Mp]]
points_B = [[Mp,Mp],[Mx-2,My-2]]
points_A = [[Mp,Mp],[Mp,My-2]]
points_H = [[Mp,Mp],[2,My-2]]

points_G = [[Mp,Mp],[2,Mp]]
points_F = [[Mp,Mp],[2,2]]
points_E = [[Mp,Mp],[Mp,2]]
points_D = [[Mp,Mp],[Mx-2,2]]

transects=[points_A,points_B,points_C,points_D,points_E,points_F,points_G,points_H]
point_names=['A','B','C','D','E','F','G','H']
dp=np.float(dkm)/np.float(resolution)
trans = ph.get_troughs(pism_infile,transects,dp)



################################################################################################
# average along profiles
profiles = {}

xav=np.zeros([len(trans),300])
yav=np.zeros([len(trans),300])
sav=np.zeros([len(trans),300])

Hav=np.zeros([Mt,len(trans),300])
mav=np.zeros([Mt,len(trans),300])
uav=np.zeros([Mt,len(trans),300])
vav=np.zeros([Mt,len(trans),300])

exp_xm = np.zeros([Mx,My])
exp_ym = np.zeros([Mx,My])

for i in range(Mx):
    exp_xm[i,:]=exp_x[i]
for j in range(My):
    exp_ym[:,j]=exp_y[j]

print('Interpolate data along profiles...')
for ti in range(Mt):
  for l,po in enumerate(trans):
    profile=[]
    for k,p in enumerate(po):
        
        i=int(np.floor(p[1]))
        j=int(np.floor(p[0]))
        di=p[0]-np.floor(p[0])
        dj=p[1]-np.floor(p[1])

        if ti==0:
            xav[l,k]=ph.interpolate_along_transect(exp_ym,i,j,di,dj)
            yav[l,k]=ph.interpolate_along_transect(exp_xm,i,j,di,dj) 
            sav[l,k]=np.sqrt(xav[l,k]**2+yav[l,k]**2)
            pair = [xav[l,k],yav[l,k],sav[l,k]]
            profile.append(pair)
        
        Hav[ti,l,k]=ph.interpolate_along_transect(exp_thk[ti],i,j,di,dj)
        uav[ti,l,k]=ph.interpolate_along_transect(exp_xvm[ti],i,j,di,dj)
        vav[ti,l,k]=ph.interpolate_along_transect(exp_yvm[ti],i,j,di,dj)
        
        #nearest neighbors
        mav[ti,l,k]=ph.nearest_along_transect(exp_mask[ti],i,j,di,dj)

    
    if ti==0:
        profiles[point_names[l]]=profile
    
########################################################################################################

if True:
  print("Write data to netCDF file "+exp_outfile)
  wrtfile = nc.Dataset(exp_outfile, 'w', format='NETCDF4_CLASSIC')

  wrtfile.createDimension('X', size=len(exp_x))
  wrtfile.createDimension('Y', size=len(exp_y))

  wrtfile.createDimension('Time1', size=Mt)
  wrtfile.createDimension('Time100', size=len(exp_t[::100]))


  nct     = wrtfile.createVariable('Time1', 'f8', ('Time1',))
  nct2    = wrtfile.createVariable('Time100', 'f8', ('Time100',))
  ncx     = wrtfile.createVariable('X', 'f8', ('X',))
  ncy     = wrtfile.createVariable('Y', 'f8', ('Y',))
    
  ncxvm   = wrtfile.createVariable('xvelmean', 'f8', ('Time100','Y', 'X'),fill_value=np.nan)
  ncyvm   = wrtfile.createVariable('yvelmean', 'f8', ('Time100','Y', 'X'),fill_value=np.nan)
  ncthk   = wrtfile.createVariable('lithk', 'f8', ('Time100','Y', 'X'),fill_value=np.nan)
  ncmask  = wrtfile.createVariable('mask', 'f8', ('Time100','Y', 'X'))
  nccrate = wrtfile.createVariable('calverate', 'f8', ('Time100','Y', 'X'),fill_value=np.nan)
  nctopg  = wrtfile.createVariable('topg', 'f8', ('Time100','Y', 'X'))
    
  ncafl  = wrtfile.createVariable('iareafl', 'f8', ('Time1'))
  ncagr  = wrtfile.createVariable('iareagr', 'f8', ('Time1'))
  nclim  = wrtfile.createVariable('lim', 'f8', ('Time1'))
  nclimw  = wrtfile.createVariable('limnsw', 'f8', ('Time1'))
  nctcf  = wrtfile.createVariable('tendlicalvf', 'f8', ('Time1'))
  nctgf  = wrtfile.createVariable('tendligroundf', 'f8', ('Time1'))
    
  wrtfile.createDimension('Profile A', size=np.shape(trans[0])[0])
  ncxa     = wrtfile.createVariable('Profile A', 'f8', ('Profile A',))
  ncsa    = wrtfile.createVariable('sA', 'f8', ('Profile A'))
  ncthka  = wrtfile.createVariable('lithkA', 'f8', ('Time1','Profile A'),fill_value=np.nan)
  ncxvma  = wrtfile.createVariable('xvelmeanA', 'f8', ('Time1','Profile A'),fill_value=np.nan)
  ncyvma  = wrtfile.createVariable('yvelmeanA', 'f8', ('Time1','Profile A'),fill_value=np.nan)
  ncmaska  = wrtfile.createVariable('maskA', 'f8', ('Time1','Profile A'))
    
  wrtfile.createDimension('Profile B', size=np.shape(trans[1])[0])
  ncxb     = wrtfile.createVariable('Profile B', 'f8', ('Profile B',))
  ncsb    = wrtfile.createVariable('sB', 'f8', ('Profile B'))
  ncthkb  = wrtfile.createVariable('lithkB', 'f8', ('Time1','Profile B'),fill_value=np.nan)
  ncxvmb  = wrtfile.createVariable('xvelmeanB', 'f8', ('Time1','Profile B'),fill_value=np.nan)
  ncyvmb  = wrtfile.createVariable('yvelmeanB', 'f8', ('Time1','Profile B'),fill_value=np.nan)
  ncmaskb  = wrtfile.createVariable('maskB', 'f8', ('Time1','Profile B'))

  wrtfile.createDimension('Profile C', size=np.shape(trans[2])[0])
  ncxc     = wrtfile.createVariable('Profile C', 'f8', ('Profile C',))
  ncsc    = wrtfile.createVariable('sC', 'f8', ('Profile C'))
  ncthkc  = wrtfile.createVariable('lithkC', 'f8', ('Time1','Profile C'),fill_value=np.nan)
  ncxvmc  = wrtfile.createVariable('xvelmeanC', 'f8', ('Time1','Profile C'),fill_value=np.nan)
  ncyvmc  = wrtfile.createVariable('yvelmeanC', 'f8', ('Time1','Profile C'),fill_value=np.nan)
  ncmaskc  = wrtfile.createVariable('maskC', 'f8', ('Time1','Profile C'))
    
  wrtfile.createDimension('Profile D', size=np.shape(trans[3])[0])
  ncxd     = wrtfile.createVariable('Profile D', 'f8', ('Profile D',))
  ncsd    = wrtfile.createVariable('sD', 'f8', ('Profile D'))
  ncthkd  = wrtfile.createVariable('lithkD', 'f8', ('Time1','Profile D'),fill_value=np.nan)
  ncxvmd  = wrtfile.createVariable('xvelmeanD', 'f8', ('Time1','Profile D'),fill_value=np.nan)
  ncyvmd  = wrtfile.createVariable('yvelmeanD', 'f8', ('Time1','Profile D'),fill_value=np.nan)
  ncmaskd  = wrtfile.createVariable('maskD', 'f8', ('Time1','Profile D'))

  wrtfile.createDimension('Profile E', size=np.shape(trans[4])[0])
  ncxe     = wrtfile.createVariable('Profile E', 'f8', ('Profile E',))
  ncse    = wrtfile.createVariable('sE', 'f8', ('Profile E'))
  ncthke  = wrtfile.createVariable('lithkE', 'f8', ('Time1','Profile E'),fill_value=np.nan)
  ncxvme  = wrtfile.createVariable('xvelmeanE', 'f8', ('Time1','Profile E'),fill_value=np.nan)
  ncyvme  = wrtfile.createVariable('yvelmeanE', 'f8', ('Time1','Profile E'),fill_value=np.nan)
  ncmaske  = wrtfile.createVariable('maskE', 'f8', ('Time1','Profile E'))
    
  wrtfile.createDimension('Profile F', size=np.shape(trans[5])[0])
  ncxf     = wrtfile.createVariable('Profile F', 'f8', ('Profile F',))
  ncsf    = wrtfile.createVariable('sF', 'f8', ('Profile F'))
  ncthkf  = wrtfile.createVariable('lithkF', 'f8', ('Time1','Profile F'),fill_value=np.nan)
  ncxvmf  = wrtfile.createVariable('xvelmeanF', 'f8', ('Time1','Profile F'),fill_value=np.nan)
  ncyvmf  = wrtfile.createVariable('yvelmeanF', 'f8', ('Time1','Profile F'),fill_value=np.nan)
  ncmaskf  = wrtfile.createVariable('maskF', 'f8', ('Time1','Profile F'))

  wrtfile.createDimension('Profile G', size=np.shape(trans[6])[0])
  ncxg     = wrtfile.createVariable('Profile G', 'f8', ('Profile G',))
  ncsg    = wrtfile.createVariable('sG', 'f8', ('Profile G'))
  ncthkg  = wrtfile.createVariable('lithkG', 'f8', ('Time1','Profile G'),fill_value=np.nan)
  ncxvmg  = wrtfile.createVariable('xvelmeanG', 'f8', ('Time1','Profile G'),fill_value=np.nan)
  ncyvmg  = wrtfile.createVariable('yvelmeanG', 'f8', ('Time1','Profile G'),fill_value=np.nan)
  ncmaskg  = wrtfile.createVariable('maskG', 'f8', ('Time1','Profile G'))
    
  wrtfile.createDimension('Profile H', size=np.shape(trans[7])[0])
  ncxh     = wrtfile.createVariable('Profile H', 'f8', ('Profile H',))
  ncsh    = wrtfile.createVariable('sH', 'f8', ('Profile H'))
  ncthkh  = wrtfile.createVariable('lithkH', 'f8', ('Time1','Profile H'),fill_value=np.nan)
  ncxvmh  = wrtfile.createVariable('xvelmeanH', 'f8', ('Time1','Profile H'),fill_value=np.nan)
  ncyvmh  = wrtfile.createVariable('yvelmeanH', 'f8', ('Time1','Profile H'),fill_value=np.nan)
  ncmaskh  = wrtfile.createVariable('maskH', 'f8', ('Time1','Profile H'))
    
  #####################################################################################
    
  nct[:] = exp_t[:]
  nct2[:]= exp_t[::100]
  ncx[:] = exp_x[:]
  ncy[:] = exp_y[:]
    

  ncxvm[:]  = exp_xvm[::100]
  ncyvm[:]  = exp_yvm[::100] 
  ncthk[:]  = exp_thk[::100]
  ncmask[:] = exp_mask[::100]
  nccrate[:]= exp_crate[::100]
  nctopg[:] = exp_topg[::100]
    

  ncafl[1:]  = exp_ts_afl[:]
  ncagr[1:]  = exp_ts_agr[:]
  nclim[1:]  = exp_ts_lim[:]
  nclimw[1:] = exp_ts_limnsw[:]
  nctcf[1:] = exp_ts_tendcf[:]
  nctgf[1:] = exp_ts_tendgf[:]
    

  cuta = np.shape(trans[0])[0]
  ncxa[:] = sav[0][0:cuta] 
  ncsa[:] = sav[0][0:cuta]
  ncthka[:] = Hav[:,0,0:cuta]
  ncmaska[:] = mav[:,0,0:cuta]
  ncxvma[:] = uav[:,0,0:cuta]
  ncyvma[:] = vav[:,0,0:cuta]
    
  cutb = np.shape(trans[1])[0]
  ncxb[:] = sav[1][0:cutb] 
  ncsb[:] = sav[1][0:cutb]
  ncthkb[:] = Hav[:,1,0:cutb]
  ncmaskb[:] = mav[:,1,0:cutb]
  ncxvmb[:] = uav[:,1,0:cutb]
  ncyvmb[:] = vav[:,1,0:cutb]
    
  cutc = np.shape(trans[2])[0]
  ncxc[:] = sav[2][0:cutc] 
  ncsc[:] = sav[2][0:cutc]
  ncthkc[:] = Hav[:,2,0:cutc]
  ncmaskc[:] = mav[:,2,0:cutc]
  ncxvmc[:] = uav[:,2,0:cutc]
  ncyvmc[:] = vav[:,2,0:cutc]

  cutd = np.shape(trans[3])[0]
  ncxd[:] = sav[3][0:cutd] 
  ncsd[:] = sav[3][0:cutd]
  ncthkd[:] = Hav[:,3,0:cutd]
  ncmaskd[:] = mav[:,3,0:cutd]
  ncxvmd[:] = uav[:,3,0:cutd]
  ncyvmd[:] = vav[:,3,0:cutd]
    
  cute = np.shape(trans[4])[0]
  ncxe[:] = sav[4][0:cute] 
  ncse[:] = sav[4][0:cute]
  ncthke[:] = Hav[:,4,0:cute]
  ncmaske[:] = mav[:,4,0:cute]
  ncxvme[:] = uav[:,4,0:cute]
  ncyvme[:] = vav[:,4,0:cute]

  cutf = np.shape(trans[5])[0]
  ncxf[:] = sav[5][0:cutf] 
  ncsf[:] = sav[5][0:cutf]
  ncthkf[:] = Hav[:,5,0:cutf]
  ncmaskf[:] = mav[:,5,0:cutf]
  ncxvmf[:] = uav[:,5,0:cutf]
  ncyvmf[:] = vav[:,5,0:cutf]
    
  cutg = np.shape(trans[6])[0]
  ncxg[:] = sav[6][0:cutg] 
  ncsg[:] = sav[6][0:cutg]
  ncthkg[:] = Hav[:,6,0:cutg]
  ncmaskg[:] = mav[:,6,0:cutg]
  ncxvmg[:] = uav[:,6,0:cutg]
  ncyvmg[:] = vav[:,6,0:cutg]

  cuth = np.shape(trans[7])[0]
  ncxh[:] = sav[7][0:cuth] 
  ncsh[:] = sav[7][0:cuth]
  ncthkh[:] = Hav[:,7,0:cuth]
  ncmaskh[:] = mav[:,7,0:cuth]
  ncxvmh[:] = uav[:,7,0:cuth]
  ncyvmh[:] = vav[:,7,0:cuth]

 
  #####################################################################################

  nct.units = 'a'
  nct2.units = 'a'
  ncx.units = 'm'
  ncy.units = 'm'
    
  ncxvm.units = 'm/a'
  ncyvm.units = 'm/a'  
  ncthk.units = 'm'
  nctopg.units = 'm'
  nccrate.units = 'm/a'
    
  ncafl.units = 'm^2'
  ncagr.units = 'm^2'
  nclim.units = 'kg'
  nclimw.units = 'kg'
  nctcf.units = 'kg/a'
  nctgf.units = 'kg/a'
    
  ncxa.units = 'm'
  ncsa.units = 'm'
  ncthka.units = 'm'
  ncxvma.units = 'm/a'
  ncyvma.units = 'm/a'
    
  ncxb.units = 'm'
  ncsb.units = 'm'
  ncthkb.units = 'm'
  ncxvmb.units = 'm/a'
  ncyvmb.units = 'm/a'
    
  ncxc.units = 'm'
  ncsc.units = 'm'
  ncthkc.units = 'm'
  ncxvmc.units = 'm/a'
  ncyvmc.units = 'm/a'
    
  ncxd.units = 'm'
  ncsd.units = 'm'
  ncthkd.units = 'm'
  ncxvmd.units = 'm/a'
  ncyvmd.units = 'm/a'
    
  ncxe.units = 'm'
  ncse.units = 'm'
  ncthke.units = 'm'
  ncxvme.units = 'm/a'
  ncyvme.units = 'm/a'
    
  ncxf.units = 'm'
  ncsf.units = 'm'
  ncthkf.units = 'm'
  ncxvmf.units = 'm/a'
  ncyvmf.units = 'm/a'
    
  ncxg.units = 'm'
  ncsg.units = 'm'
  ncthkg.units = 'm'
  ncxvmg.units = 'm/a'
  ncyvmg.units = 'm/a'
    
  ncxh.units = 'm'
  ncsh.units = 'm'
  ncthkh.units = 'm'
  ncxvmh.units = 'm/a'
  ncyvmh.units = 'm/a'

    
  #####################################################################################


  ncxvm.Standard_name   = 'land_ice_vertical_mean_x_velocity'
  ncyvm.Standard_name   = 'land_ice_vertical_mean_y_velocity'
  ncthk.Standard_name   = 'land_ice_thickness'
  nctopg.Standard_name  = 'bedrock_altitude'
  nccrate.Standard_name = 'calving_rate'
    
  ncafl.Standard_name  = 'floating_ice_shelf_area'
  ncagr.Standard_name  = 'grounded_ice_sheet_area'
  nclim.Standard_name  =  'land_ice_mass'
  nclimw.Standard_name = 'land_ice_mass_not_displacing_sea_water'
  nctcf.Standard_name  =  'tendency_of_land_ice_mass_due_to_calving'
  nctgf.Standard_name  =  'tendency_of_grounded_ice_mass'
    
  ncsa.Standard_name   = 'distance_along_profile_A'
  ncthka.Standard_name = 'land_ice_thickness_along_profile_A'
  ncxvma.Standard_name = 'land_ice_vertical_mean_x_velocity_along_profile_A'
  ncyvma.Standard_name = 'land_ice_vertical_mean_y_velocity_along_profile_A'
    
  ncsb.Standard_name   = 'distance_along_profile_B'
  ncthkb.Standard_name = 'land_ice_thickness_along_profile_B'
  ncxvmb.Standard_name = 'land_ice_vertical_mean_x_velocity_along_profile_B'
  ncyvmb.Standard_name = 'land_ice_vertical_mean_y_velocity_along_profile_B'    
    
  ncsc.Standard_name   = 'distance_along_profile_C'
  ncthkc.Standard_name = 'land_ice_thickness_along_profile_C'
  ncxvmc.Standard_name = 'land_ice_vertical_mean_x_velocity_along_profile_C'
  ncyvmc.Standard_name = 'land_ice_vertical_mean_y_velocity_along_profile_C'
    
  ncsd.Standard_name   = 'distance_along_profile_D'
  ncthkd.Standard_name = 'land_ice_thickness_along_profile_D'
  ncxvmd.Standard_name = 'land_ice_vertical_mean_x_velocity_along_profile_D'
  ncyvmd.Standard_name = 'land_ice_vertical_mean_y_velocity_along_profile_D'
    
  ncse.Standard_name   = 'distance_along_profile_E'
  ncthke.Standard_name = 'land_ice_thickness_along_profile_E'
  ncxvme.Standard_name = 'land_ice_vertical_mean_x_velocity_along_profile_E'
  ncyvme.Standard_name = 'land_ice_vertical_mean_y_velocity_along_profile_E'
    
  ncsf.Standard_name   = 'distance_along_profile_F'
  ncthkf.Standard_name = 'land_ice_thickness_along_profile_F'
  ncxvmf.Standard_name = 'land_ice_vertical_mean_x_velocity_along_profile_F'
  ncyvmf.Standard_name = 'land_ice_vertical_mean_y_velocity_along_profile_F'
    
  ncsg.Standard_name   = 'distance_along_profile_G'
  ncthkg.Standard_name = 'land_ice_thickness_along_profile_G'
  ncxvmg.Standard_name = 'land_ice_vertical_mean_x_velocity_along_profile_G'
  ncyvmg.Standard_name = 'land_ice_vertical_mean_y_velocity_along_profile_G'
    
  ncsh.Standard_name   = 'distance_along_profile_H'
  ncthkh.Standard_name = 'land_ice_thickness_along_profile_H'
  ncxvmh.Standard_name = 'land_ice_vertical_mean_x_velocity_along_profile_H'
  ncyvmh.Standard_name = 'land_ice_vertical_mean_y_velocity_along_profile_H'
    
  #####################################################################################
    
  ncmask.flag_values = '1, 2, 3'
  ncmask.flag_meanings = '1=grounded ice, 2=floating ice, 3=open ocean'
   
  ncmaska.flag_values = '1, 2, 3'
  ncmaska.flag_meanings = '1=grounded ice, 2=floating ice, 3=open ocean'

  ncmaskb.flag_values = '1, 2, 3'
  ncmaskb.flag_meanings = '1=grounded ice, 2=floating ice, 3=open ocean'

  ncmaskc.flag_values = '1, 2, 3'
  ncmaskc.flag_meanings = '1=grounded ice, 2=floating ice, 3=open ocean'

  ncmaskd.flag_values = '1, 2, 3'
  ncmaskd.flag_meanings = '1=grounded ice, 2=floating ice, 3=open ocean'

  ncmaske.flag_values = '1, 2, 3'
  ncmaske.flag_meanings = '1=grounded ice, 2=floating ice, 3=open ocean'

  ncmaskf.flag_values = '1, 2, 3'
  ncmaskf.flag_meanings = '1=grounded ice, 2=floating ice, 3=open ocean'

  ncmaskg.flag_values = '1, 2, 3'
  ncmaskg.flag_meanings = '1=grounded ice, 2=floating ice, 3=open ocean'

  ncmaskh.flag_values = '1, 2, 3'
  ncmaskh.flag_meanings = '1=grounded ice, 2=floating ice, 3=open ocean'
    

  now = datetime.datetime.now().strftime("%B %d, %Y")
  wrtfile.comment  = "CalvingMIP contribution postprocessed by torsten.albrecht@pik-potsdam.de at " + now
  wrtfile.institution = 'Potsdam Institute for Climate Impact Research (PIK), Germany'
  wrtfile.inputdata = 'PISM code from https://github.com/pism/pism/tree/pik/calving_rate_given'
    
  wrtfile.close()


  print('Data successfully saved to', exp_outfile)

