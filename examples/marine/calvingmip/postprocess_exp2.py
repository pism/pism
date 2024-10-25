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
vcr=1e-8 # terminal velocity threshold


pismpath     = "/p/tmp/albrecht/pism23/calvmip/circular/exp2-05km-dir-fr0-retreat/"
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
exp_subh = datnc.variables["ice_area_specific_volume"][:]
datnc.close()

exp_thks=exp_thk+exp_subh # sum of thk and ice_area_specific_volume
exp_thkm=np.ma.masked_array(exp_thk,mask=exp_thk<=0.0) # exclude empty cells from mean
exp_subhm=np.ma.masked_array(exp_subh,mask=exp_subh<=0.0) # exclude empty subgrid cells

exp_xvmc=np.ma.masked_array(exp_xvm,mask=np.abs(exp_xvm)<=vcr) # remove small values
exp_yvmc=np.ma.masked_array(exp_yvm,mask=np.abs(exp_yvm)<=vcr)


#print(exp_t)
Mt,Mx,My = np.shape(exp_mask)

exp_thk[exp_thk == 0.0] = np.nan
exp_xvm[exp_xvm == 0.0] = np.nan
exp_yvm[exp_yvm == 0.0] = np.nan

# avoid empty first snapshot
exp_xvm[0]=exp_xvm[1]
exp_yvm[0]=exp_yvm[1]

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
#print(Mx,My,Mp)

points_C = [[Mp,Mp],[Mx-2,Mp]]
points_B = [[Mp,Mp],[Mx-2,My-2]]
points_A = [[Mp,Mp],[Mp,My-2]]
points_H = [[Mp,Mp],[1,My-2]]

points_G = [[Mp,Mp],[1,Mp]]
points_F = [[Mp,Mp],[1,1]]
points_E = [[Mp,Mp],[Mp,1]]
points_D = [[Mp,Mp],[Mx-2,1]]

transects=[points_A,points_B,points_C,points_D,points_E,points_F,points_G,points_H]
point_names=['A','B','C','D','E','F','G','H']
dp=np.float(dkm)/np.float(resolution)
trans = ph.get_troughs(pism_infile,transects,dp)



################################################################################################
# average along profiles
profiles = {}

lentr=len(trans)

xav=np.zeros([lentr,300])
yav=np.zeros([lentr,300])
sav=np.zeros([lentr,300])

Hav=np.zeros([Mt,lentr,300])
mav=np.zeros([Mt,lentr,300])
uav=np.zeros([Mt,lentr,300])
vav=np.zeros([Mt,lentr,300])

Hcf=np.zeros([Mt,lentr])
xcf=np.zeros([Mt,lentr])
ycf=np.zeros([Mt,lentr])
vxcf=np.zeros([Mt,lentr])
vycf=np.zeros([Mt,lentr])

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
    dx=po[1][1]-po[0][1] # from 5km interpolation along transects, defining a direction
    dy=po[1][0]-po[0][0]

    for k,p in enumerate(po):
        
        i=int(np.floor(p[1]))
        j=int(np.floor(p[0]))
        dj=p[0]-np.floor(p[0])
        di=p[1]-np.floor(p[1])

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


        #re-define for calving front interpolation
        i=int(np.around(p[1]))
        j=int(np.around(p[0]))

        if (exp_mask[ti,int(np.sign(dx)+i),int(np.sign(dy)+j)]==3.0 and exp_mask[ti,i,j]==2.0):

          Atot      = 9.0 # area of neighbor cells around i,j
          Hcf[ti,l] = np.nanmean(exp_thkm[ti,i-1:i+2,j-1:j+2]) # mean over cells with thk>0
          Acov      = np.sum((exp_thks[ti,i-1:i+2,j-1:j+2])/Hcf[ti,l]) # area covered by cells with thk>0 or icea_area_specific_volume>0
          ds        = 3.0*(Acov/Atot-0.5) # average distance of i,j from ice front, assuming half of 9 cells would be filled with thk>0 
          ycf[ti,l] = exp_ym[i,j]+ds*dy*(exp_ym[i,j+1]-exp_ym[i,j]) # direction is determined by dy
          xcf[ti,l] = exp_xm[i,j]+ds*dx*(exp_xm[i+1,j]-exp_xm[i,j]) # direction is determined vy dx

          vxcf[ti,l]= np.mean(exp_xvmc[ti,i-1:i+2,j-1:j+2])
          vycf[ti,l]= np.mean(exp_yvmc[ti,i-1:i+2,j-1:j+2])


    if ti==0:
        profiles[point_names[l]]=profile


################################################################################################
# get ice covered quadrants

Anw=np.zeros([Mt])
Ane=np.zeros([Mt])
Asw=np.zeros([Mt])
Ase=np.zeros([Mt])

Axy=(resolution*1e3)**2 #m2

print('Aggregate ice covered area...')

for ti,t in enumerate(exp_t):
    for i,xi in enumerate(exp_x):
      for j,yj in enumerate(exp_y):
        if exp_mask[ti,i,j]<3.0:

            if xi>=0 and yj>=0:
                if xi==0.0 and yj==0.0:
                    Ane[ti]+=0.25
                elif xi==0.0 or yj==0.0:
                    Ane[ti]+=0.5
                else:
                    Ane[ti]+=1.0
            if xi>=0 and yj<=0:
                if xi==0.0 and yj==0.0:
                    Anw[ti]+=0.25
                elif xi==0.0 or yj==0.0:
                    Anw[ti]+=0.5
                else:
                    Anw[ti]+=1.0
            if xi<=0 and yj>=0:
                if xi==0.0 and yj==0.0:
                    Ase[ti]+=0.25
                elif xi==0.0 or yj==0.0:
                    Ase[ti]+=0.5
                else:
                    Ase[ti]+=1.0
            if xi<=0 and yj<=0:
                if xi==0.0 and yj==0.0:
                    Asw[ti]+=0.25
                elif xi==0.0 or yj==0.0:
                    Asw[ti]+=0.5
                else:
                    Asw[ti]+=1.0
                
#print(Ane,Anw,Ase,Asw)

Anw*=Axy
Ane*=Axy
Asw*=Axy
Ase*=Axy


###########################################################################################################

lent=Mt

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


  ncatnw  = wrtfile.createVariable('iareatotalNW', 'f8', ('Time1'))
  ncatne  = wrtfile.createVariable('iareatotalNE', 'f8', ('Time1'))
  ncatsw  = wrtfile.createVariable('iareatotalSW', 'f8', ('Time1'))
  ncatse  = wrtfile.createVariable('iareatotalSE', 'f8', ('Time1'))

  ncxcfa  = wrtfile.createVariable('xcfA', 'f8', ('Time1'),fill_value=np.nan)
  ncycfa  = wrtfile.createVariable('ycfA', 'f8', ('Time1'),fill_value=np.nan)
  ncxvmcfa  = wrtfile.createVariable('xvelmeancfA', 'f8', ('Time1'))
  ncyvmcfa  = wrtfile.createVariable('yvelmeancfA', 'f8', ('Time1'))
  ncthkcfa  = wrtfile.createVariable('lithkcfA', 'f8', ('Time1'))

  ncxcfb  = wrtfile.createVariable('xcfB', 'f8', ('Time1'),fill_value=np.nan)
  ncycfb  = wrtfile.createVariable('ycfB', 'f8', ('Time1'),fill_value=np.nan)
  ncxvmcfb  = wrtfile.createVariable('xvelmeancfB', 'f8', ('Time1'))
  ncyvmcfb  = wrtfile.createVariable('yvelmeancfB', 'f8', ('Time1'))
  ncthkcfb  = wrtfile.createVariable('lithkcfB', 'f8', ('Time1'))

  ncxcfc  = wrtfile.createVariable('xcfC', 'f8', ('Time1'),fill_value=np.nan)
  ncycfc  = wrtfile.createVariable('ycfC', 'f8', ('Time1'),fill_value=np.nan)
  ncxvmcfc  = wrtfile.createVariable('xvelmeancfC', 'f8', ('Time1'))
  ncyvmcfc  = wrtfile.createVariable('yvelmeancfC', 'f8', ('Time1'))
  ncthkcfc  = wrtfile.createVariable('lithkcfC', 'f8', ('Time1'))

  ncxcfd  = wrtfile.createVariable('xcfD', 'f8', ('Time1'),fill_value=np.nan)
  ncycfd  = wrtfile.createVariable('ycfD', 'f8', ('Time1'),fill_value=np.nan)
  ncxvmcfd  = wrtfile.createVariable('xvelmeancfD', 'f8', ('Time1'))
  ncyvmcfd  = wrtfile.createVariable('yvelmeancfD', 'f8', ('Time1'))
  ncthkcfd  = wrtfile.createVariable('lithkcfD', 'f8', ('Time1'))

  ncxcfe  = wrtfile.createVariable('xcfE', 'f8', ('Time1'), fill_value=np.nan)
  ncycfe  = wrtfile.createVariable('ycfE', 'f8', ('Time1'),fill_value=np.nan)
  ncxvmcfe  = wrtfile.createVariable('xvelmeancfE', 'f8', ('Time1'))
  ncyvmcfe  = wrtfile.createVariable('yvelmeancfE', 'f8', ('Time1'))
  ncthkcfe  = wrtfile.createVariable('lithkcfE', 'f8', ('Time1'))
 
  ncxcff  = wrtfile.createVariable('xcfF', 'f8', ('Time1'),fill_value=np.nan)
  ncycff  = wrtfile.createVariable('ycfF', 'f8', ('Time1'),fill_value=np.nan)
  ncxvmcff  = wrtfile.createVariable('xvelmeancfF', 'f8', ('Time1'))
  ncyvmcff  = wrtfile.createVariable('yvelmeancfF', 'f8', ('Time1'))
  ncthkcff  = wrtfile.createVariable('lithkcfF', 'f8', ('Time1'))

  ncxcfg  = wrtfile.createVariable('xcfG', 'f8', ('Time1'),fill_value=np.nan)
  ncycfg  = wrtfile.createVariable('ycfG', 'f8', ('Time1'),fill_value=np.nan)
  ncxvmcfg  = wrtfile.createVariable('xvelmeancfG', 'f8', ('Time1'))
  ncyvmcfg  = wrtfile.createVariable('yvelmeancfG', 'f8', ('Time1'))
  ncthkcfg  = wrtfile.createVariable('lithkcfG', 'f8', ('Time1'))

  ncxcfh  = wrtfile.createVariable('xcfH', 'f8', ('Time1'),fill_value=np.nan)
  ncycfh  = wrtfile.createVariable('ycfH', 'f8', ('Time1'),fill_value=np.nan)
  ncxvmcfh  = wrtfile.createVariable('xvelmeancfH', 'f8', ('Time1'))
  ncyvmcfh  = wrtfile.createVariable('yvelmeancfH', 'f8', ('Time1'))
  ncthkcfh  = wrtfile.createVariable('lithkcfH', 'f8', ('Time1'))

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


  ncatnw[:] = Anw[:]
  ncatne[:] = Ane[:]
  ncatsw[:] = Asw[:]
  ncatse[:] = Ase[:]

  ncxcfa[:] = xcf[:lent,0]
  ncycfa[:] = ycf[:lent,0]
  ncxvmcfa[:] = vxcf[:lent,0]
  ncyvmcfa[:] = vycf[:lent,0]
  ncthkcfa[:] = Hcf[:lent,0]

  ncxcfb[:] = xcf[:lent,1]
  ncycfb[:] = ycf[:lent,1]
  ncxvmcfb[:] = vxcf[:lent,1]
  ncyvmcfb[:] = vycf[:lent,1]
  ncthkcfb[:] = Hcf[:lent,1]

  ncxcfc[:] = xcf[:lent,2]
  ncycfc[:] = ycf[:lent,2]
  ncxvmcfc[:] = vxcf[:lent,2]
  ncyvmcfc[:] = vycf[:lent,2]
  ncthkcfc[:] = Hcf[:lent,2]

  ncxcfd[:] = xcf[:lent,3]
  ncycfd[:] = ycf[:lent,3]
  ncxvmcfd[:] = vxcf[:lent,3]
  ncyvmcfd[:] = vycf[:lent,3]
  ncthkcfd[:] = Hcf[:lent,3]

  ncxcfe[:] = xcf[:lent,4]
  ncycfe[:] = ycf[:lent,4]
  ncxvmcfe[:] = vxcf[:lent,4]
  ncyvmcfe[:] = vycf[:lent,4]
  ncthkcfe[:] = Hcf[:lent,4]

  ncxcff[:] = xcf[:lent,5]
  ncycff[:] = ycf[:lent,5]
  ncxvmcff[:] = vxcf[:lent,5]
  ncyvmcff[:] = vycf[:lent,5]
  ncthkcff[:] = Hcf[:lent,5]

  ncxcfg[:] = xcf[:lent,6]
  ncycfg[:] = ycf[:lent,6]
  ncxvmcfg[:] = vxcf[:lent,6]
  ncyvmcfg[:] = vycf[:lent,6]
  ncthkcfg[:] = Hcf[:lent,6]

  ncxcfh[:] = xcf[:lent,7]
  ncycfh[:] = ycf[:lent,7]
  ncxvmcfh[:] = vxcf[:lent,7]
  ncyvmcfh[:] = vycf[:lent,7]
  ncthkcfh[:] = Hcf[:lent,7]


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


  ncatnw.units = 'm^2'
  ncatne.units = 'm^2'
  ncatsw.units = 'm^2'
  ncatse.units = 'm^2'

  ncxcfa.units = 'm'
  ncycfa.units = 'm'
  ncxvmcfa.units = 'm/a'
  ncyvmcfa.units = 'm/a'
  ncthkcfa.units = 'm'

  ncxcfb.units = 'm'
  ncycfb.units = 'm'
  ncxvmcfb.units = 'm/a'
  ncyvmcfb.units = 'm/a'
  ncthkcfb.units = 'm'
 
  ncxcfc.units = 'm'
  ncycfc.units = 'm'
  ncxvmcfc.units = 'm/a'
  ncyvmcfc.units = 'm/a'
  ncthkcfc.units = 'm'

  ncxcfd.units = 'm'
  ncycfd.units = 'm'
  ncxvmcfd.units = 'm/a'
  ncyvmcfd.units = 'm/a'
  ncthkcfd.units = 'm'

  ncxcfe.units = 'm'
  ncycfe.units = 'm'
  ncxvmcfe.units = 'm/a'
  ncyvmcfe.units = 'm/a'
  ncthkcfe.units = 'm'

  ncxcff.units = 'm'
  ncycff.units = 'm'
  ncxvmcff.units = 'm/a'
  ncyvmcff.units = 'm/a'
  ncthkcff.units = 'm'

  ncxcfg.units = 'm'
  ncycfg.units = 'm'
  ncxvmcfg.units = 'm/a'
  ncyvmcfg.units = 'm/a'
  ncthkcfg.units = 'm'

  ncxcfh.units = 'm'
  ncycfh.units = 'm'
  ncxvmcfh.units = 'm/a'
  ncyvmcfh.units = 'm/a'
  ncthkcfh.units = 'm'


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


  ncatnw.Standard_name  = 'total_ice_area_NorthWest'
  ncatne.Standard_name  = 'total_ice_area_NorthEast'
  ncatsw.Standard_name  = 'total_ice_area_SouthWest'
  ncatse.Standard_name  = 'total_ice_area_SouthEast'

  ncxcfa.Standard_name   = 'x_calving_front_on_profile_A'
  ncycfa.Standard_name   = 'y_calving_front_on_profile_A'
  ncxvmcfa.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_profile_A'
  ncyvmcfa.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_profile_A'
  ncthkcfa.Standard_name = 'land_ice_thickness_at_calving_front_on_profile_A'

  ncxcfb.Standard_name   = 'x_calving_front_on_profile_B'
  ncycfb.Standard_name   = 'y_calving_front_on_profile_B'
  ncxvmcfb.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_profile_B'
  ncyvmcfb.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_profile_B'
  ncthkcfb.Standard_name = 'land_ice_thickness_at_calving_front_on_profile_B'

  ncxcfc.Standard_name   = 'x_calving_front_on_profile_C'
  ncycfc.Standard_name   = 'y_calving_front_on_profile_C'
  ncxvmcfc.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_profile_C'
  ncyvmcfc.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_profile_C'
  ncthkcfc.Standard_name = 'land_ice_thickness_at_calving_front_on_profile_C'

  ncxcfd.Standard_name   = 'x_calving_front_on_profile_D'
  ncycfd.Standard_name   = 'y_calving_front_on_profile_D'
  ncxvmcfd.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_profile_D'
  ncyvmcfd.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_profile_D'
  ncthkcfd.Standard_name = 'land_ice_thickness_at_calving_front_on_profile_D'

  ncxcfe.Standard_name   = 'x_calving_front_on_profile_E'
  ncycfe.Standard_name   = 'y_calving_front_on_profile_e'
  ncxvmcfe.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_profile_E'
  ncyvmcfe.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_profile_E'
  ncthkcfe.Standard_name = 'land_ice_thickness_at_calving_front_on_profile_E'

  ncxcff.Standard_name   = 'x_calving_front_on_profile_F'
  ncycff.Standard_name   = 'y_calving_front_on_profile_F'
  ncxvmcff.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_profile_F'
  ncyvmcff.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_profile_F'
  ncthkcff.Standard_name = 'land_ice_thickness_at_calving_front_on_profile_F'

  ncxcfg.Standard_name   = 'x_calving_front_on_profile_G'
  ncycfg.Standard_name   = 'y_calving_front_on_profile_G'
  ncxvmcfg.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_profile_G'
  ncyvmcfg.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_profile_G'
  ncthkcfg.Standard_name = 'land_ice_thickness_at_calving_front_on_profile_G'

  ncxcfh.Standard_name   = 'x_calving_front_on_profile_H'
  ncycfh.Standard_name   = 'y_calving_front_on_profile_H'
  ncxvmcfh.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_profile_H'
  ncyvmcfh.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_profile_H'
  ncthkcfh.Standard_name = 'land_ice_thickness_at_calving_front_on_profile_H'

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

