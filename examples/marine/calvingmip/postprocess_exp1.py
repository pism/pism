#!/usr/bin/env python3

# Creates output file from PISM EXP1 result to upload for CalvinMIP, 
# as instructed from https://github.com/JRowanJordan/CalvingMIP/wiki/Experiment-1

import numpy as np
import netCDF4 as nc

secperyear=365*24*3600
rhoi_bm = 917.0
rhosw_bm = 1028.0

resolution=5.0
dkm=5.0 #km steps


pism_outfile = "/p/tmp/albrecht/pism23/calvmip/circular/exp1-05km-dir/results/extra_exp1_test.nc" 
pism_tsfile  = "/p/tmp/albrecht/pism23/calvmip/circular/exp1-05km-dir/results/ts_exp1_last.nc"
pism_infile = "/p/tmp/albrecht/pism23/calvmip/circular/exp1-05km-dir/input/circular_input_5km.nc"

exp1_outfile = "CalvingMIP_EXP1_PISM_PIK.nc"
#postprocessfile = "/p/tmp/albrecht/pism23/calvmip/circular/exp1-05km-dir/postprocess/CalvingMIP_EXP1_PISM_PIK.nc"

#################################################################################################

datnc = nc.Dataset(pism_outfile,"r")
exp1_x = datnc.variables["x"][:]
exp1_y = datnc.variables["y"][:]
exp1_t = [datnc.variables["time"][-1]/secperyear-110000.0]
#exp1_tb = datnc.variables["time_bounds"][:]//secperyear-110000.0


exp1_xvm = datnc.variables["xvelmean"][-1,:]*secperyear
exp1_yvm = datnc.variables["yvelmean"][-1,:]*secperyear
exp1_thk = datnc.variables["lithk"][-1,:]
exp1_mask = datnc.variables["mask"][-1,:]-1
exp1_crate = datnc.variables["calvingmip_calving_rate"][-1:]*secperyear
exp1_topg = datnc.variables["topg"][-1,:]

#exp1_afl = datnc.variables["sftflf"][:]
#exp1_agr = datnc.variables["sftgrf"][:]

datnc.close()

Mx,My = np.shape(exp1_mask)




datnc = nc.Dataset(pism_tsfile,"r")

exp1_ts_t = datnc.variables["time"][:]/secperyear-110000.0
#exp1_ts_tb = datnc.variables["time_bounds"][:]/secperyear-110000.0

exp1_ts_afl = datnc.variables["iareafl"][:]
exp1_ts_agr = datnc.variables["iareagr"][:]
exp1_ts_lim = datnc.variables["lim"][:]
exp1_ts_limnsw = datnc.variables["limnsw"][:]
exp1_ts_tendcf = datnc.variables["tendlicalvf"][:]*secperyear
exp1_ts_tendgf = datnc.variables["tendligroundf"][:]*secperyear

datnc.close()

###############################################################################################

# define piecewise linear transects between i-j-locations 
def get_troughs(of,trans):

  try:
    compfile = nc.Dataset(of, 'r')
  except Exception:
    print("Specify NetCDF input file.")
    print(of)
    exit(-1)
  xobs = np.squeeze(compfile.variables["x"][:])
  yobs = np.squeeze(compfile.variables["y"][:])
  #mobs = np.squeeze(compfile.variables["mask"][:])
  compfile.close()

  trans_points=[]
  
  for l,points in enumerate(trans):

    pkm=[]
    for pp in points:
      pkm.append([xobs[pp[0]],yobs[pp[1]]])

    dp=np.float(dkm)/np.float(resolution)
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
trans = get_troughs(pism_infile,transects)


# plot profiles
if False:
    fig = plt.figure(4,figsize=(8,8))
    ax = fig.add_subplot(1,1,1)

    for cnt,p in enumerate(transects):
        p1=p[0]
        p2=p[1]
    
        print(p,np.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2))
        ax.plot([p1[0],p2[0]],[p1[1],p2[1]],label=point_names[cnt])
      
    ax.set_xlim([0,Mx])
    ax.set_ylim([0,My])
    ax.legend()


# average along profiles
profiles = {}

xav=np.zeros([len(trans),300])
yav=np.zeros([len(trans),300])

Hav=np.zeros([len(trans),300])
mav=np.zeros([len(trans),300])
uav=np.zeros([len(trans),300])
vav=np.zeros([len(trans),300])
sav=np.zeros([len(trans),300])

exp1_xm = np.zeros([Mx,My])
exp1_ym = np.zeros([Mx,My])

for i in range(Mx):
    exp1_xm[i,:]=exp1_x[i]
for j in range(My):
    exp1_ym[:,j]=exp1_y[j]
    

for l,po in enumerate(trans):
    profile=[]
    for k,p in enumerate(po):
        
        i=int(np.floor(p[1]))
        j=int(np.floor(p[0]))
        di=p[0]-np.floor(p[0])
        dj=p[1]-np.floor(p[1])
        
        xav[l,k]=interpolate_along_transect(exp1_ym,i,j,di,dj)
        yav[l,k]=interpolate_along_transect(exp1_xm,i,j,di,dj) 
        sav[l,k]=np.sqrt(xav[l,k]**2+yav[l,k]**2)
        
        Hav[l,k]=interpolate_along_transect(exp1_thk,i,j,di,dj)
        mav[l,k]=interpolate_along_transect(exp1_mask,i,j,di,dj)
        uav[l,k]=interpolate_along_transect(exp1_xvm,i,j,di,dj)
        vav[l,k]=interpolate_along_transect(exp1_yvm,i,j,di,dj)
     
        pair = [xav[l,k],yav[l,k],sav[l,k]]
        profile.append(pair)
    
    profiles[point_names[l]]=profile
    
########################################################################################################

if True:
  print("Write to netcdf file "+exp1_outfile)
  wrtfile = nc.Dataset(exp1_outfile, 'w', format='NETCDF4_CLASSIC')
  wrtfile.createDimension('X', size=len(exp1_x))
  wrtfile.createDimension('Y', size=len(exp1_y))
  #wrtfile.createDimension('Time', size=None)
  wrtfile.createDimension('Time', size=1)
  wrtfile.createDimension('nv', size=2)

  nct     = wrtfile.createVariable('Time', 'f8', ('Time',))
  #nctb    = wrtfile.createVariable('Time_bnds', 'f8', ('Time','nv'))
  ncx     = wrtfile.createVariable('X', 'f8', ('X',))
  ncy     = wrtfile.createVariable('Y', 'f8', ('Y',))
    
  ncxvm   = wrtfile.createVariable('xvelmean', 'f8', ('Time','Y', 'X'),fill_value=np.nan)
  ncyvm   = wrtfile.createVariable('yvelmean', 'f8', ('Time','Y', 'X'),fill_value=np.nan)
  ncthk   = wrtfile.createVariable('lithk', 'f8', ('Time','Y', 'X'),fill_value=np.nan)
  ncmask  = wrtfile.createVariable('mask', 'f8', ('Time','Y', 'X'))
  nccrate = wrtfile.createVariable('calverate', 'f8', ('Time','Y', 'X'),fill_value=np.nan)
  nctopg  = wrtfile.createVariable('topg', 'f8', ('Time','Y', 'X'))
    
  ncafl  = wrtfile.createVariable('iareafl', 'f8', ('Time'))
  ncagr  = wrtfile.createVariable('iareagr', 'f8', ('Time'))
  nclim  = wrtfile.createVariable('lim', 'f8', ('Time'))
  nclimw  = wrtfile.createVariable('limnsw', 'f8', ('Time'))
  nctcf  = wrtfile.createVariable('tendlicalvf', 'f8', ('Time'))
  nctgf  = wrtfile.createVariable('tendligroundf', 'f8', ('Time'))
    
  exp1_thk[exp1_thk == 0.0] = np.nan
  exp1_xvm[exp1_xvm == 0.0] = np.nan
  exp1_yvm[exp1_yvm == 0.0] = np.nan

  wrtfile.createDimension('Profile A', size=np.shape(trans[0])[0])
  ncxa     = wrtfile.createVariable('Profile A', 'f8', ('Profile A',))
  ncsa    = wrtfile.createVariable('sA', 'f8', ('Profile A'))
  ncthka  = wrtfile.createVariable('lithkA', 'f8', ('Profile A'),fill_value=np.nan)
  ncxvma  = wrtfile.createVariable('xvelmeanA', 'f8', ('Profile A'),fill_value=np.nan)
  ncyvma  = wrtfile.createVariable('yvelmeanA', 'f8', ('Profile A'),fill_value=np.nan)
  ncmaska  = wrtfile.createVariable('maskA', 'f8', ('Profile A'))
    
  wrtfile.createDimension('Profile B', size=np.shape(trans[1])[0])
  ncxb     = wrtfile.createVariable('Profile B', 'f8', ('Profile B',))
  ncsb    = wrtfile.createVariable('sB', 'f8', ('Profile B'))
  ncthkb  = wrtfile.createVariable('lithkB', 'f8', ('Profile B'),fill_value=np.nan)
  ncxvmb  = wrtfile.createVariable('xvelmeanB', 'f8', ('Profile B'),fill_value=np.nan)
  ncyvmb  = wrtfile.createVariable('yvelmeanB', 'f8', ('Profile B'),fill_value=np.nan)
  ncmaskb  = wrtfile.createVariable('maskB', 'f8', ('Profile B'))

  wrtfile.createDimension('Profile C', size=np.shape(trans[2])[0])
  ncxc     = wrtfile.createVariable('Profile C', 'f8', ('Profile C',))
  ncsc    = wrtfile.createVariable('sC', 'f8', ('Profile C'))
  ncthkc  = wrtfile.createVariable('lithkC', 'f8', ('Profile C'),fill_value=np.nan)
  ncxvmc  = wrtfile.createVariable('xvelmeanC', 'f8', ('Profile C'),fill_value=np.nan)
  ncyvmc  = wrtfile.createVariable('yvelmeanC', 'f8', ('Profile C'),fill_value=np.nan)
  ncmaskc  = wrtfile.createVariable('maskC', 'f8', ('Profile C'))
    
  wrtfile.createDimension('Profile D', size=np.shape(trans[3])[0])
  ncxd     = wrtfile.createVariable('Profile D', 'f8', ('Profile D',))
  ncsd    = wrtfile.createVariable('sD', 'f8', ('Profile D'))
  ncthkd  = wrtfile.createVariable('lithkD', 'f8', ('Profile D'),fill_value=np.nan)
  ncxvmd  = wrtfile.createVariable('xvelmeanD', 'f8', ('Profile D'),fill_value=np.nan)
  ncyvmd  = wrtfile.createVariable('yvelmeanD', 'f8', ('Profile D'),fill_value=np.nan)
  ncmaskd  = wrtfile.createVariable('maskD', 'f8', ('Profile D'))

  wrtfile.createDimension('Profile E', size=np.shape(trans[4])[0])
  ncxe     = wrtfile.createVariable('Profile E', 'f8', ('Profile E',))
  ncse    = wrtfile.createVariable('sE', 'f8', ('Profile E'))
  ncthke  = wrtfile.createVariable('lithkE', 'f8', ('Profile E'),fill_value=np.nan)
  ncxvme  = wrtfile.createVariable('xvelmeanE', 'f8', ('Profile E'),fill_value=np.nan)
  ncyvme  = wrtfile.createVariable('yvelmeanE', 'f8', ('Profile E'),fill_value=np.nan)
  ncmaske  = wrtfile.createVariable('maskE', 'f8', ('Profile E'))
    
  wrtfile.createDimension('Profile F', size=np.shape(trans[5])[0])
  ncxf     = wrtfile.createVariable('Profile F', 'f8', ('Profile F',))
  ncsf    = wrtfile.createVariable('sF', 'f8', ('Profile F'))
  ncthkf  = wrtfile.createVariable('lithkF', 'f8', ('Profile F'),fill_value=np.nan)
  ncxvmf  = wrtfile.createVariable('xvelmeanF', 'f8', ('Profile F'),fill_value=np.nan)
  ncyvmf  = wrtfile.createVariable('yvelmeanF', 'f8', ('Profile F'),fill_value=np.nan)
  ncmaskf  = wrtfile.createVariable('maskF', 'f8', ('Profile F'))

  wrtfile.createDimension('Profile G', size=np.shape(trans[6])[0])
  ncxg     = wrtfile.createVariable('Profile G', 'f8', ('Profile G',))
  ncsg    = wrtfile.createVariable('sG', 'f8', ('Profile G'))
  ncthkg  = wrtfile.createVariable('lithkG', 'f8', ('Profile G'),fill_value=np.nan)
  ncxvmg  = wrtfile.createVariable('xvelmeanG', 'f8', ('Profile G'),fill_value=np.nan)
  ncyvmg  = wrtfile.createVariable('yvelmeanG', 'f8', ('Profile G'),fill_value=np.nan)
  ncmaskg  = wrtfile.createVariable('maskG', 'f8', ('Profile G'))
    
  wrtfile.createDimension('Profile H', size=np.shape(trans[7])[0])
  ncxh     = wrtfile.createVariable('Profile H', 'f8', ('Profile H',))
  ncsh    = wrtfile.createVariable('sH', 'f8', ('Profile H'))
  ncthkh  = wrtfile.createVariable('lithkH', 'f8', ('Profile H'),fill_value=np.nan)
  ncxvmh  = wrtfile.createVariable('xvelmeanH', 'f8', ('Profile H'),fill_value=np.nan)
  ncyvmh  = wrtfile.createVariable('yvelmeanH', 'f8', ('Profile H'),fill_value=np.nan)
  ncmaskh  = wrtfile.createVariable('maskH', 'f8', ('Profile H'))
    
  #####################################################################################
    
  nct[:] = exp1_t[:]
  ncx[:] = exp1_x[:]
  ncy[:] = exp1_y[:]
  #nctb[:] = exp1_tb[:]
    
  ncxvm[:]  = exp1_xvm[:]
  ncyvm[:]  = exp1_yvm[:] 
  ncthk[:]  = exp1_thk[:]
  ncmask[:] = exp1_mask[:] #mask		Ice mask	grounded=1, floating=2, open ocean=3
  nccrate[:]= exp1_crate[:]
  nctopg[:] = exp1_topg[:]
    
  ncafl[:]  = exp1_ts_afl[:]
  ncagr[:]  = exp1_ts_agr[:]
  nclim[:]  = exp1_ts_lim[:]
  nclimw[:] = exp1_ts_limnsw[:]
  nctcf[:] = exp1_ts_tendcf[:]
  nctgf[:] = exp1_ts_tendgf[:]
    
  cuta = np.shape(trans[0])[0]
  ncxa[:] = sav[0][0:cuta] 
  ncsa[:] = sav[0][0:cuta]
  ncthka[:] = Hav[0][0:cuta]
  ncmaska[:] = mav[0][0:cuta]
  ncxvma[:] = uav[0][0:cuta]
  ncyvma[:] = vav[0][0:cuta]
    
  cutb = np.shape(trans[1])[0]
  ncxb[:] = sav[1][0:cutb] 
  ncsb[:] = sav[1][0:cutb]
  ncthkb[:] = Hav[1][0:cutb]
  ncmaskb[:] = mav[1][0:cutb]
  ncxvmb[:] = uav[1][0:cutb]
  ncyvmb[:] = vav[1][0:cutb]
    
  cutc = np.shape(trans[2])[0]
  ncxc[:] = sav[2][0:cutc] 
  ncsc[:] = sav[2][0:cutc]
  ncthkc[:] = Hav[2][0:cutc]
  ncmaskc[:] = mav[2][0:cutc]
  ncxvmc[:] = uav[2][0:cutc]
  ncyvmc[:] = vav[2][0:cutc]

  cutd = np.shape(trans[3])[0]
  ncxd[:] = sav[3][0:cutd] 
  ncsd[:] = sav[3][0:cutd]
  ncthkd[:] = Hav[3][0:cutd]
  ncmaskd[:] = mav[3][0:cutd]
  ncxvmd[:] = uav[3][0:cutd]
  ncyvmd[:] = vav[3][0:cutd]
    
  cute = np.shape(trans[4])[0]
  ncxe[:] = sav[4][0:cute] 
  ncse[:] = sav[4][0:cute]
  ncthke[:] = Hav[4][0:cute]
  ncmaske[:] = mav[4][0:cute]
  ncxvme[:] = uav[4][0:cute]
  ncyvme[:] = vav[4][0:cute]

  cutf = np.shape(trans[5])[0]
  ncxf[:] = sav[5][0:cutf] 
  ncsf[:] = sav[5][0:cutf]
  ncthkf[:] = Hav[5][0:cutf]
  ncmaskf[:] = mav[5][0:cutf]
  ncxvmf[:] = uav[5][0:cutf]
  ncyvmf[:] = vav[5][0:cutf]
    
  cutg = np.shape(trans[6])[0]
  ncxg[:] = sav[6][0:cutg] 
  ncsg[:] = sav[6][0:cutg]
  ncthkg[:] = Hav[6][0:cutg]
  ncmaskg[:] = mav[6][0:cutg]
  ncxvmg[:] = uav[6][0:cutg]
  ncyvmg[:] = vav[6][0:cutg]

  cuth = np.shape(trans[7])[0]
  ncxh[:] = sav[7][0:cuth] 
  ncsh[:] = sav[7][0:cuth]
  ncthkh[:] = Hav[7][0:cuth]
  ncmaskh[:] = mav[7][0:cuth]
  ncxvmh[:] = uav[7][0:cuth]
  ncyvmh[:] = vav[7][0:cuth]

 
  #####################################################################################

  nct.units = 'a'
  #nctb.units = 'a'
  ncx.units = 'm'
  ncy.units = 'm'
    
  ncxvm.units = 'm/a'
  ncyvm.units = 'm/a'  
  ncthk.units = 'm'
  nctopg.units = 'm'
  #ncmask.units = ''
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
  #ncmaska.units = ''
  ncxvma.units = 'm/a'
  ncyvma.units = 'm/a'
    
  ncxb.units = 'm'
  ncsb.units = 'm'
  ncthkb.units = 'm'
  #ncmaskb.units = ''
  ncxvmb.units = 'm/a'
  ncyvmb.units = 'm/a'
    
  ncxc.units = 'm'
  ncsc.units = 'm'
  ncthkc.units = 'm'
  #ncmaskc.units = ''
  ncxvmc.units = 'm/a'
  ncyvmc.units = 'm/a'
    
  ncxd.units = 'm'
  ncsd.units = 'm'
  ncthkd.units = 'm'
  #ncmaskd.units = ''
  ncxvmd.units = 'm/a'
  ncyvmd.units = 'm/a'
    
  ncxe.units = 'm'
  ncse.units = 'm'
  ncthke.units = 'm'
  #ncmaske.units = ''
  ncxvme.units = 'm/a'
  ncyvme.units = 'm/a'
    
  ncxf.units = 'm'
  ncsf.units = 'm'
  ncthkf.units = 'm'
  #ncmaskf.units = ''
  ncxvmf.units = 'm/a'
  ncyvmf.units = 'm/a'
    
  ncxg.units = 'm'
  ncsg.units = 'm'
  ncthkg.units = 'm'
  #ncmaskg.units = ''
  ncxvmg.units = 'm/a'
  ncyvmg.units = 'm/a'
    
  ncxh.units = 'm'
  ncsh.units = 'm'
  ncthkh.units = 'm'
  #ncmaskh.units = ''
  ncxvmh.units = 'm/a'
  ncyvmh.units = 'm/a'

    
  #####################################################################################

  #nct.Standard_name = 'a'
  #nctb.Standard_name = 'a'
  #ncx.Standard_name = "projection_x_coordinate"
  #ncy.Standard_name = "projection_y_coordinate"
    
  ncxvm.Standard_name   = 'land_ice_vertical_mean_x_velocity'
  ncyvm.Standard_name   = 'land_ice_vertical_mean_y_velocity'
  ncthk.Standard_name   = 'land_ice_thickness'
  nctopg.Standard_name  = 'bedrock_altitude'
  #ncmask.Standard_name = ''
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
    

  import datetime
  now = datetime.datetime.now().strftime("%B %d, %Y")
  wrtfile.comment  = "CalvingMIP contribution postprocessed by torsten.albrecht@pik-potsdam.de at " + now
  wrtfile.institution = 'Potsdam Institute for Climate Impact Research (PIK), Germany'
  wrtfile.inputdata = 'PISM code from https://github.com/pism/pism/tree/pik/calving_rate_given'
    
  wrtfile.close()

  print('Data successfully saved to', exp1_outfile)

