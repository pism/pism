#!/usr/bin/env python3

# Creates output file from PISM EXP1 result to upload for CalvinMIP, 
# as instructed from https://github.com/JRowanJordan/CalvingMIP/wiki/Experiment-1

import numpy as np
import netCDF4 as nc
import datetime
import postprocess_helper as ph

secperyear=365*24*3600

resolution=5.0
dkm=5.0 #km steps

pismpath     = "/p/tmp/albrecht/pism23/calvmip/circular/exp1-05km-dir/"
pism_outfile = pismpath + "results/result_exp1c.nc"
pism_tsfile  = pismpath + "results/ts_exp1c.nc"

pism_infile = pismpath + "input/circular_input_5km.nc"

exp_outfile = "CalvingMIP_EXP1_PISM_PIK.nc"

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
exp_mask  = datnc.variables["mask"][:]-1
exp_crate = datnc.variables["calvingmip_calving_rate"][:]*secperyear
exp_topg  = datnc.variables["topg"][:]
exp_subh = datnc.variables["ice_area_specific_volume"][:]
datnc.close()

#print(exp_t)
Mt,Mx,My = np.shape(exp_mask)

exp_thk[exp_thk == 0.0] = np.nan
exp_xvm[exp_xvm == 0.0] = np.nan
exp_yvm[exp_yvm == 0.0] = np.nan

datnc = nc.Dataset(pism_tsfile,"r")
exp_ts_t      = datnc.variables["time"][:]/secperyear-120000.0
exp_ts_afl    = datnc.variables["iareafl"][:]
exp_ts_agr    = datnc.variables["iareagr"][:]
exp_ts_lim    = datnc.variables["lim"][:]
exp_ts_limnsw = datnc.variables["limnsw"][:]
exp_ts_tendcf = datnc.variables["tendlicalvf"][:]*secperyear
exp_ts_tendgf = datnc.variables["tendligroundf"][:]*secperyear
datnc.close()
#print(exp_ts_t)

###############################################################################################

# define profiles
Mp=int((Mx-1)/2.0)
#print(Mx,My,Mp)

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

Hcf =np.zeros([Mt,len(trans)])
xcf =np.zeros([Mt,len(trans)])
ycf =np.zeros([Mt,len(trans)])
vxcf=np.zeros([Mt,len(trans)])
vycf=np.zeros([Mt,len(trans)])

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


        # get marginal values along profiles
        if l==0:
            if exp_mask[ti,i+1,j]==3.0 and exp_mask[ti,i,j]==2.0:
              Hcf[ti,l]=exp_thk[ti,i,j]
              vxcf[ti,l]=exp_xvm[ti,i,j]
              vycf[ti,l]=exp_yvm[ti,i,j]
              ycf[ti,l]=exp_ym[i,j]+(0.5+exp_subh[ti,i+1,j]/Hcf[ti,l])*(exp_ym[i+1,j]-exp_ym[i,j])
              xcf[ti,l]=exp_xm[i,j]+(0.5+exp_subh[ti,i+1,j]/Hcf[ti,l])*(exp_xm[i+1,j]-exp_xm[i,j])
        elif l==2:
            if exp_mask[ti,i,j+1]==3.0 and exp_mask[ti,i,j]==2.0:
              Hcf[ti,l]=exp_thk[ti,i,j]
              vxcf[ti,l]=exp_xvm[ti,i,j]
              vycf[ti,l]=exp_yvm[ti,i,j]
              ycf[ti,l]=exp_ym[i,j]+(0.5+exp_subh[ti,i,j+1]/Hcf[ti,l])*(exp_ym[i,j+1]-exp_ym[i,j])
              xcf[ti,l]=exp_xm[i,j]+(0.5+exp_subh[ti,i,j+1]/Hcf[ti,l])*(exp_xm[i,j+1]-exp_xm[i,j])
        elif l==4:
            if exp_mask[ti,i-1,j]==3.0 and exp_mask[ti,i,j]==2.0:
              Hcf[ti,l]=exp_thk[ti,i,j]
              vxcf[ti,l]=exp_xvm[ti,i,j]
              vycf[ti,l]=exp_yvm[ti,i,j]
              ycf[ti,l]=exp_ym[i,j]+(0.5+exp_subh[ti,i-1,j]/Hcf[ti,l])*(exp_ym[i-1,j]-exp_ym[i,j])
              xcf[ti,l]=exp_xm[i,j]+(0.5+exp_subh[ti,i-1,j]/Hcf[ti,l])*(exp_xm[i-1,j]-exp_xm[i,j])
        elif l==6:
            if (exp_mask[ti,i,j-1]==3.0 and exp_mask[ti,i,j]==2.0):
              Hcf[ti,l]=exp_thk[ti,i,j]
              vxcf[ti,l]=exp_xvm[ti,i,j]
              vycf[ti,l]=exp_yvm[ti,i,j]
              ycf[ti,l]=exp_ym[i,j]+(0.5+exp_subh[ti,i,j-1]/Hcf[ti,l])*(exp_ym[i,j-1]-exp_ym[i,j])
              xcf[ti,l]=exp_xm[i,j]+(0.5+exp_subh[ti,i,j-1]/Hcf[ti,l])*(exp_xm[i,j-1]-exp_xm[i,j])
        elif l==1:
            if (exp_mask[ti,i-1,j+1]==3.0 and exp_mask[ti,i,j]==2.0):
              Hcf[ti,l]=exp_thk[ti,i,j]
              vxcf[ti,l]=exp_xvm[ti,i,j]
              vycf[ti,l]=exp_yvm[ti,i,j]
              if exp_mask[ti,i,j+1]==3.0:
                ds=0.5
                if exp_subh[ti,i,j+1]<=Hcf[ti,l]/2.0:
                  ds+=np.sqrt(exp_subh[ti,i,j+1]/(2.0*Hcf[ti,l])) #rectangular triangle
                elif exp_subh[ti,i,j+1]>Hcf[ti,l]/2.0:
                  ds+=np.sqrt(exp_subh[ti,i,j+1]/(2.0*Hcf[ti,l])-0.25)+0.5
              else: #exp_mask[ti,i,j+1]==2.0
                ds=1.5
                if exp_subh[ti,i-1,j+1]<=Hcf[ti,l]/2.0:
                  ds+=np.sqrt(exp_subh[ti,i-1,j+1]/(2.0*Hcf[ti,l])) #rectangular triangle
                elif exp_subh[ti,i-1,j+1]>Hcf[ti,l]/2.0:
                  ds+=np.sqrt(exp_subh[ti,i-1,j+1]/(2.0*Hcf[ti,l])-0.25)+0.5
              ycf[ti,l]=exp_ym[i,j]+0.5*ds*(exp_ym[i-1,j+1]-exp_ym[i,j])
              xcf[ti,l]=exp_xm[i,j]+0.5*ds*(exp_xm[i-1,j+1]-exp_xm[i,j])

        elif l==3:
            if (exp_mask[ti,i-1,j+1]==3.0 and exp_mask[ti,i,j]==2.0):
              Hcf[ti,l]=exp_thk[ti,i,j]
              vxcf[ti,l]=exp_xvm[ti,i,j]
              vycf[ti,l]=exp_yvm[ti,i,j]
              if exp_mask[ti,i,j+1]==3.0:
                ds=0.5
                if exp_subh[ti,i,j+1]<=Hcf[ti,l]/2.0:
                  ds+=np.sqrt(exp_subh[ti,i,j+1]/(2.0*Hcf[ti,l])) #rectangular triangle
                elif exp_subh[ti,i,j+1]>Hcf[ti,l]/2.0:
                  ds+=np.sqrt(exp_subh[ti,i,j+1]/(2.0*Hcf[ti,l])-0.25)+0.5
              else: #exp_mask[ti,i,j+1]==2.0
                ds=1.5
                if exp_subh[ti,i-1,j+1]<=Hcf[ti,l]/2.0:
                  ds+=np.sqrt(exp_subh[ti,i-1,j+1]/(2.0*Hcf[ti,l])) #rectangular triangle
                elif exp_subh[ti,i-1,j+1]>Hcf[ti,l]/2.0:
                  ds+=np.sqrt(exp_subh[ti,i-1,j+1]/(2.0*Hcf[ti,l])-0.25)+0.5
              ycf[ti,l]=exp_ym[i,j]+0.5*ds*(exp_ym[i-1,j+1]-exp_ym[i,j])
              xcf[ti,l]=exp_xm[i,j]+0.5*ds*(exp_xm[i-1,j+1]-exp_xm[i,j])

        elif l==5:
            if (exp_mask[ti,i-1,j-1]==3.0 and exp_mask[ti,i,j]==2.0):
              Hcf[ti,l]=exp_thk[ti,i,j]
              vxcf[ti,l]=exp_xvm[ti,i,j]
              vycf[ti,l]=exp_yvm[ti,i,j]
              if exp_mask[ti,i,j-1]==3.0:
                ds=0.5
                if exp_subh[ti,i,j-1]<=Hcf[ti,l]/2.0:
                  ds+=np.sqrt(exp_subh[ti,i,j-1]/(2.0*Hcf[ti,l])) #rectangular triangle
                elif exp_subh[ti,i,j-1]>Hcf[ti,l]/2.0:
                  ds+=np.sqrt(exp_subh[ti,i,j-1]/(2.0*Hcf[ti,l])-0.25)+0.5
              else: #exp_mask[ti,i,j-1]==2.0
                ds=1.5
                if exp_subh[ti,i-1,j-1]<=Hcf[ti,l]/2.0:
                  ds+=np.sqrt(exp_subh[ti,i-1,j-1]/(2.0*Hcf[ti,l])) #rectangular triangle
                elif exp_subh[ti,i-1,j-1]>Hcf[ti,l]/2.0:
                  ds+=np.sqrt(exp_subh[ti,i-1,j-1]/(2.0*Hcf[ti,l])-0.25)+0.5
              ycf[ti,l]=exp_ym[i,j]+0.5*ds*(exp_ym[i-1,j-1]-exp_ym[i,j])
              xcf[ti,l]=exp_xm[i,j]+0.5*ds*(exp_xm[i-1,j-1]-exp_xm[i,j])
 
        elif l==7:
            if (exp_mask[ti,i+1,j-1]==3.0 and exp_mask[ti,i,j]==2.0):
              Hcf[ti,l]=exp_thk[ti,i,j]
              vxcf[ti,l]=exp_xvm[ti,i,j]
              vycf[ti,l]=exp_yvm[ti,i,j]
              if exp_mask[ti,i,j-1]==3.0:
                ds=0.5
                if exp_subh[ti,i,j-1]<=Hcf[ti,l]/2.0:
                  ds+=np.sqrt(exp_subh[ti,i,j-1]/(2.0*Hcf[ti,l])) #rectangular triangle
                elif exp_subh[ti,i,j-1]>Hcf[ti,l]/2.0:
                  ds+=np.sqrt(exp_subh[ti,i,j-1]/(2.0*Hcf[ti,l])-0.25)+0.5
              else: #exp_mask[ti,i,j-1]==2.0
                ds=1.5
                if exp_subh[ti,i+1,j-1]<=Hcf[ti,l]/2.0:
                  ds+=np.sqrt(exp_subh[ti,i+1,j-1]/(2.0*Hcf[ti,l])) #rectangular triangle
                elif exp_subh[ti,i+1,j-1]>Hcf[ti,l]/2.0:
                  ds+=np.sqrt(exp_subh[ti,i+1,j-1]/(2.0*Hcf[ti,l])-0.25)+0.5
              ycf[ti,l]=exp_ym[i,j]+0.5*ds*(exp_ym[i+1,j-1]-exp_ym[i,j])
              xcf[ti,l]=exp_xm[i,j]+0.5*ds*(exp_xm[i+1,j-1]-exp_xm[i,j])

    if ti==0:
        profiles[point_names[l]]=profile


################################################################################################
# get ice covered quadrants

Anw=0.0 #+1-1
Ane=0.0 #+1+1
Asw=0.0 #-1-1
Ase=0.0 #-1+1

Axy=(resolution*1e3)**2 #m2

for i,xi in enumerate(exp_x):
    for j,yj in enumerate(exp_y):
        if exp_mask[0,i,j]<3.0:
            if xi>=0 and yj>=0:
                # 0,0 is in the middle of the center grid cell
                if xi==0.0 and yj==0.0:
                    Ane+=0.25
                elif xi==0.0 or yj==0.0:
                    Ane+=0.5
                else:
                    Ane+=1.0
            if xi>=0 and yj<=0:
                if xi==0.0 and yj==0.0:
                    Anw+=0.25
                elif xi==0.0 or yj==0.0:
                    Anw+=0.5
                else:
                    Anw+=1.0
            if xi<=0 and yj>=0:
                if xi==0.0 and yj==0.0:
                    Ase+=0.25
                elif xi==0.0 or yj==0.0:
                    Ase+=0.5
                else:
                    Ase+=1.0
            if xi<=0 and yj<=0:
                if xi==0.0 and yj==0.0:
                    Asw+=0.25
                elif xi==0.0 or yj==0.0:
                    Asw+=0.5
                else:
                    Asw+=1.0


Anw*=Axy
Ane*=Axy
Asw*=Axy
Ase*=Axy

print(Ane,Anw,Ase,Asw)


########################################################################################################

if True:
  print("Write data to netCDF file "+exp_outfile)
  wrtfile = nc.Dataset(exp_outfile, 'w', format='NETCDF4_CLASSIC')

  wrtfile.createDimension('X', size=len(exp_x))
  wrtfile.createDimension('Y', size=len(exp_y))

  wrtfile.createDimension('Time1', size=1)
  #wrtfile.createDimension('Time100', size=len(exp_t[::100]))


  nct     = wrtfile.createVariable('Time1', 'f8', ('Time1',))
  #nct2    = wrtfile.createVariable('Time100', 'f8', ('Time100',))
  ncx     = wrtfile.createVariable('X', 'f8', ('X',))
  ncy     = wrtfile.createVariable('Y', 'f8', ('Y',))

  ncxvm   = wrtfile.createVariable('xvelmean', 'f8', ('Time1','Y', 'X'),fill_value=np.nan)
  ncyvm   = wrtfile.createVariable('yvelmean', 'f8', ('Time1','Y', 'X'),fill_value=np.nan)
  ncthk   = wrtfile.createVariable('lithk', 'f8', ('Time1','Y', 'X'),fill_value=np.nan)
  ncmask  = wrtfile.createVariable('mask', 'f8', ('Time1','Y', 'X'))
  nccrate = wrtfile.createVariable('calverate', 'f8', ('Time1','Y', 'X'),fill_value=np.nan)
  nctopg  = wrtfile.createVariable('topg', 'f8', ('Time1','Y', 'X'))

  ncafl  = wrtfile.createVariable('iareafl', 'f8', ('Time1'))
  ncagr  = wrtfile.createVariable('iareagr', 'f8', ('Time1'))
  nclim  = wrtfile.createVariable('lim', 'f8', ('Time1'))
  nclimw  = wrtfile.createVariable('limnsw', 'f8', ('Time1'))
  nctcf  = wrtfile.createVariable('tendlicalvf', 'f8', ('Time1'))
  nctgf  = wrtfile.createVariable('tendligroundf', 'f8', ('Time1'))

  ncatnw  = wrtfile.createVariable('iareatotalNW', 'f8', ('Time1'))
  ncatne  = wrtfile.createVariable('iareatotalNE', 'f8', ('Time1'))
  ncatsw  = wrtfile.createVariable('iareatotalSW', 'f8', ('Time1'))
  ncatse  = wrtfile.createVariable('iareatotalSE', 'f8', ('Time1'))

  wrtfile.createDimension('Profile A', size=np.shape(trans[0])[0])
  ncxa     = wrtfile.createVariable('Profile A', 'f8', ('Profile A',))
  ncsa    = wrtfile.createVariable('sA', 'f8', ('Profile A'))
  ncthka  = wrtfile.createVariable('lithkA', 'f8', ('Time1','Profile A'),fill_value=np.nan)
  ncxvma  = wrtfile.createVariable('xvelmeanA', 'f8', ('Time1','Profile A'),fill_value=np.nan)
  ncyvma  = wrtfile.createVariable('yvelmeanA', 'f8', ('Time1','Profile A'),fill_value=np.nan)
  ncmaska  = wrtfile.createVariable('maskA', 'f8', ('Time1','Profile A'))

  ncxcfa  = wrtfile.createVariable('xcfA', 'f8',fill_value=np.nan)
  ncycfa  = wrtfile.createVariable('ycfA', 'f8',fill_value=np.nan)
  ncxvmcfa  = wrtfile.createVariable('xvelmeancfA', 'f8')
  ncyvmcfa  = wrtfile.createVariable('yvelmeancfA', 'f8')
  ncthkcfa  = wrtfile.createVariable('lithkcfA', 'f8')

  wrtfile.createDimension('Profile B', size=np.shape(trans[1])[0])
  ncxb     = wrtfile.createVariable('Profile B', 'f8', ('Profile B',))
  ncsb    = wrtfile.createVariable('sB', 'f8', ('Profile B'))
  ncthkb  = wrtfile.createVariable('lithkB', 'f8', ('Time1','Profile B'),fill_value=np.nan)
  ncxvmb  = wrtfile.createVariable('xvelmeanB', 'f8', ('Time1','Profile B'),fill_value=np.nan)
  ncyvmb  = wrtfile.createVariable('yvelmeanB', 'f8', ('Time1','Profile B'),fill_value=np.nan)
  ncmaskb  = wrtfile.createVariable('maskB', 'f8', ('Time1','Profile B'))

  ncxcfb  = wrtfile.createVariable('xcfB', 'f8',fill_value=np.nan)
  ncycfb  = wrtfile.createVariable('ycfB', 'f8',fill_value=np.nan)
  ncxvmcfb  = wrtfile.createVariable('xvelmeancfB', 'f8')
  ncyvmcfb  = wrtfile.createVariable('yvelmeancfB', 'f8')
  ncthkcfb  = wrtfile.createVariable('lithkcfB', 'f8')

  wrtfile.createDimension('Profile C', size=np.shape(trans[2])[0])
  ncxc     = wrtfile.createVariable('Profile C', 'f8', ('Profile C',))
  ncsc    = wrtfile.createVariable('sC', 'f8', ('Profile C'))
  ncthkc  = wrtfile.createVariable('lithkC', 'f8', ('Time1','Profile C'),fill_value=np.nan)
  ncxvmc  = wrtfile.createVariable('xvelmeanC', 'f8', ('Time1','Profile C'),fill_value=np.nan)
  ncyvmc  = wrtfile.createVariable('yvelmeanC', 'f8', ('Time1','Profile C'),fill_value=np.nan)
  ncmaskc  = wrtfile.createVariable('maskC', 'f8', ('Time1','Profile C'))

  ncxcfc  = wrtfile.createVariable('xcfC', 'f8',fill_value=np.nan)
  ncycfc  = wrtfile.createVariable('ycfC', 'f8',fill_value=np.nan)
  ncxvmcfc  = wrtfile.createVariable('xvelmeancfC', 'f8')
  ncyvmcfc  = wrtfile.createVariable('yvelmeancfC', 'f8')
  ncthkcfc  = wrtfile.createVariable('lithkcfC', 'f8')

  wrtfile.createDimension('Profile D', size=np.shape(trans[3])[0])
  ncxd     = wrtfile.createVariable('Profile D', 'f8', ('Profile D',))
  ncsd    = wrtfile.createVariable('sD', 'f8', ('Profile D'))
  ncthkd  = wrtfile.createVariable('lithkD', 'f8', ('Time1','Profile D'),fill_value=np.nan)
  ncxvmd  = wrtfile.createVariable('xvelmeanD', 'f8', ('Time1','Profile D'),fill_value=np.nan)
  ncyvmd  = wrtfile.createVariable('yvelmeanD', 'f8', ('Time1','Profile D'),fill_value=np.nan)
  ncmaskd  = wrtfile.createVariable('maskD', 'f8', ('Time1','Profile D'))

  ncxcfd  = wrtfile.createVariable('xcfD', 'f8',fill_value=np.nan)
  ncycfd  = wrtfile.createVariable('ycfD', 'f8',fill_value=np.nan)
  ncxvmcfd  = wrtfile.createVariable('xvelmeancfD', 'f8')
  ncyvmcfd  = wrtfile.createVariable('yvelmeancfD', 'f8')
  ncthkcfd  = wrtfile.createVariable('lithkcfD', 'f8')

  wrtfile.createDimension('Profile E', size=np.shape(trans[4])[0])
  ncxe     = wrtfile.createVariable('Profile E', 'f8', ('Profile E',))
  ncse    = wrtfile.createVariable('sE', 'f8', ('Profile E'))
  ncthke  = wrtfile.createVariable('lithkE', 'f8', ('Time1','Profile E'),fill_value=np.nan)
  ncxvme  = wrtfile.createVariable('xvelmeanE', 'f8', ('Time1','Profile E'),fill_value=np.nan)
  ncyvme  = wrtfile.createVariable('yvelmeanE', 'f8', ('Time1','Profile E'),fill_value=np.nan)
  ncmaske  = wrtfile.createVariable('maskE', 'f8', ('Time1','Profile E'))

  ncxcfe  = wrtfile.createVariable('xcfE', 'f8',fill_value=np.nan)
  ncycfe  = wrtfile.createVariable('ycfE', 'f8',fill_value=np.nan)
  ncxvmcfe  = wrtfile.createVariable('xvelmeancfE', 'f8')
  ncyvmcfe  = wrtfile.createVariable('yvelmeancfE', 'f8')
  ncthkcfe  = wrtfile.createVariable('lithkcfE', 'f8')

  wrtfile.createDimension('Profile F', size=np.shape(trans[5])[0])
  ncxf     = wrtfile.createVariable('Profile F', 'f8', ('Profile F',))
  ncsf    = wrtfile.createVariable('sF', 'f8', ('Profile F'))
  ncthkf  = wrtfile.createVariable('lithkF', 'f8', ('Time1','Profile F'),fill_value=np.nan)
  ncxvmf  = wrtfile.createVariable('xvelmeanF', 'f8', ('Time1','Profile F'),fill_value=np.nan)
  ncyvmf  = wrtfile.createVariable('yvelmeanF', 'f8', ('Time1','Profile F'),fill_value=np.nan)
  ncmaskf  = wrtfile.createVariable('maskF', 'f8', ('Time1','Profile F'))

  ncxcff  = wrtfile.createVariable('xcfF', 'f8',fill_value=np.nan)
  ncycff  = wrtfile.createVariable('ycfF', 'f8',fill_value=np.nan)
  ncxvmcff  = wrtfile.createVariable('xvelmeancfF', 'f8')
  ncyvmcff  = wrtfile.createVariable('yvelmeancfF', 'f8')
  ncthkcff  = wrtfile.createVariable('lithkcfF', 'f8')

  wrtfile.createDimension('Profile G', size=np.shape(trans[6])[0])
  ncxg     = wrtfile.createVariable('Profile G', 'f8', ('Profile G',))
  ncsg    = wrtfile.createVariable('sG', 'f8', ('Profile G'))
  ncthkg  = wrtfile.createVariable('lithkG', 'f8', ('Time1','Profile G'),fill_value=np.nan)
  ncxvmg  = wrtfile.createVariable('xvelmeanG', 'f8', ('Time1','Profile G'),fill_value=np.nan)
  ncyvmg  = wrtfile.createVariable('yvelmeanG', 'f8', ('Time1','Profile G'),fill_value=np.nan)
  ncmaskg  = wrtfile.createVariable('maskG', 'f8', ('Time1','Profile G'))

  ncxcfg  = wrtfile.createVariable('xcfG', 'f8',fill_value=np.nan)
  ncycfg  = wrtfile.createVariable('ycfG', 'f8',fill_value=np.nan)
  ncxvmcfg  = wrtfile.createVariable('xvelmeancfG', 'f8')
  ncyvmcfg  = wrtfile.createVariable('yvelmeancfG', 'f8')
  ncthkcfg  = wrtfile.createVariable('lithkcfG', 'f8')

  wrtfile.createDimension('Profile H', size=np.shape(trans[7])[0])
  ncxh     = wrtfile.createVariable('Profile H', 'f8', ('Profile H',))
  ncsh    = wrtfile.createVariable('sH', 'f8', ('Profile H'))
  ncthkh  = wrtfile.createVariable('lithkH', 'f8', ('Time1','Profile H'),fill_value=np.nan)
  ncxvmh  = wrtfile.createVariable('xvelmeanH', 'f8', ('Time1','Profile H'),fill_value=np.nan)
  ncyvmh  = wrtfile.createVariable('yvelmeanH', 'f8', ('Time1','Profile H'),fill_value=np.nan)
  ncmaskh  = wrtfile.createVariable('maskH', 'f8', ('Time1','Profile H'))

  ncxcfh  = wrtfile.createVariable('xcfH', 'f8',fill_value=np.nan)
  ncycfh  = wrtfile.createVariable('ycfH', 'f8',fill_value=np.nan)
  ncxvmcfh  = wrtfile.createVariable('xvelmeancfH', 'f8')
  ncyvmcfh  = wrtfile.createVariable('yvelmeancfH', 'f8')
  ncthkcfh  = wrtfile.createVariable('lithkcfH', 'f8')


  #####################################################################################

  nct[:] = exp_t[:]
  #nct2[:]= exp_t[::100]
  ncx[:] = exp_x[:]
  ncy[:] = exp_y[:]


  ncxvm[:]  = exp_xvm[:]
  ncyvm[:]  = exp_yvm[:] 
  ncthk[:]  = exp_thk[:]
  ncmask[:] = exp_mask[:]
  nccrate[:]= exp_crate[:]
  nctopg[:] = exp_topg[:]

  ncafl[:]  = exp_ts_afl[-1]
  ncagr[:]  = exp_ts_agr[-1]
  nclim[:]  = exp_ts_lim[-1]
  nclimw[:] = exp_ts_limnsw[-1]
  nctcf[:] = exp_ts_tendcf[-1]
  nctgf[:] = exp_ts_tendgf[-1]

  ncatnw[:] = Anw
  ncatne[:] = Ane
  ncatsw[:] = Asw
  ncatse[:] = Ase

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


  ncxcfa[:] = xcf[:,0]
  ncycfa[:] = ycf[:,0]
  ncxvmcfa[:] = vxcf[:,0]
  ncyvmcfa[:] = vycf[:,0]
  ncthkcfa[:] = Hcf[:,0]

  ncxcfb[:] = xcf[:,1]
  ncycfb[:] = ycf[:,1]
  ncxvmcfb[:] = vxcf[:,1]
  ncyvmcfb[:] = vycf[:,1]
  ncthkcfb[:] = Hcf[:,1]

  ncxcfc[:] = xcf[:,2]
  ncycfc[:] = ycf[:,2]
  ncxvmcfc[:] = vxcf[:,2]
  ncyvmcfc[:] = vycf[:,2]
  ncthkcfc[:] = Hcf[:,2]

  ncxcfd[:] = xcf[:,3]
  ncycfd[:] = ycf[:,3]
  ncxvmcfd[:] = vxcf[:,3]
  ncyvmcfd[:] = vycf[:,3]
  ncthkcfd[:] = Hcf[:,3]

  ncxcfe[:] = xcf[:,4]
  ncycfe[:] = ycf[:,4]
  ncxvmcfe[:] = vxcf[:,4]
  ncyvmcfe[:] = vycf[:,4]
  ncthkcfe[:] = Hcf[:,4]

  ncxcff[:] = xcf[:,5]
  ncycff[:] = ycf[:,5]
  ncxvmcff[:] = vxcf[:,5]
  ncyvmcff[:] = vycf[:,5]
  ncthkcff[:] = Hcf[:,5]

  ncxcfg[:] = xcf[:,6]
  ncycfg[:] = ycf[:,6]
  ncxvmcfg[:] = vxcf[:,6]
  ncyvmcfg[:] = vycf[:,6]
  ncthkcfg[:] = Hcf[:,6]

  ncxcfh[:] = xcf[:,7]
  ncycfh[:] = ycf[:,7]
  ncxvmcfh[:] = vxcf[:,7]
  ncyvmcfh[:] = vycf[:,7]
  ncthkcfh[:] = Hcf[:,7]

 
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

  ncatnw.units = 'm^2'
  ncatne.units = 'm^2'
  ncatsw.units = 'm^2'
  ncatse.units = 'm^2'

  ncxa.units = 'm'
  ncsa.units = 'm'
  ncthka.units = 'm'
  #ncmaska.units = ''
  ncxvma.units = 'm/a'
  ncyvma.units = 'm/a'

  ncxcfa.units = 'm'
  ncycfa.units = 'm'
  ncxvmcfa.units = 'm/a'
  ncyvmcfa.units = 'm/a'
  ncthkcfa.units = 'm'

  ncxb.units = 'm'
  ncsb.units = 'm'
  ncthkb.units = 'm'
  #ncmaskb.units = ''
  ncxvmb.units = 'm/a'
  ncyvmb.units = 'm/a'

  ncxcfb.units = 'm'
  ncycfb.units = 'm'
  ncxvmcfb.units = 'm/a'
  ncyvmcfb.units = 'm/a'
  ncthkcfb.units = 'm'

  ncxc.units = 'm'
  ncsc.units = 'm'
  ncthkc.units = 'm'
  #ncmaskc.units = ''
  ncxvmc.units = 'm/a'
  ncyvmc.units = 'm/a'

  ncxcfc.units = 'm'
  ncycfc.units = 'm'
  ncxvmcfc.units = 'm/a'
  ncyvmcfc.units = 'm/a'
  ncthkcfc.units = 'm'

  ncxd.units = 'm'
  ncsd.units = 'm'
  ncthkd.units = 'm'
  #ncmaskd.units = ''
  ncxvmd.units = 'm/a'
  ncyvmd.units = 'm/a'

  ncxcfd.units = 'm'
  ncycfd.units = 'm'
  ncxvmcfd.units = 'm/a'
  ncyvmcfd.units = 'm/a'
  ncthkcfd.units = 'm'

  ncxe.units = 'm'
  ncse.units = 'm'
  ncthke.units = 'm'
  #ncmaske.units = ''
  ncxvme.units = 'm/a'
  ncyvme.units = 'm/a'

  ncxcfe.units = 'm'
  ncycfe.units = 'm'
  ncxvmcfe.units = 'm/a'
  ncyvmcfe.units = 'm/a'
  ncthkcfe.units = 'm'

  ncxf.units = 'm'
  ncsf.units = 'm'
  ncthkf.units = 'm'
  #ncmaskf.units = ''
  ncxvmf.units = 'm/a'
  ncyvmf.units = 'm/a'

  ncxcff.units = 'm'
  ncycff.units = 'm'
  ncxvmcff.units = 'm/a'
  ncyvmcff.units = 'm/a'
  ncthkcff.units = 'm'

  ncxg.units = 'm'
  ncsg.units = 'm'
  ncthkg.units = 'm'
  #ncmaskg.units = ''
  ncxvmg.units = 'm/a'
  ncyvmg.units = 'm/a'

  ncxcfg.units = 'm'
  ncycfg.units = 'm'
  ncxvmcfg.units = 'm/a'
  ncyvmcfg.units = 'm/a'
  ncthkcfg.units = 'm'

  ncxh.units = 'm'
  ncsh.units = 'm'
  ncthkh.units = 'm'
  #ncmaskh.units = ''
  ncxvmh.units = 'm/a'
  ncyvmh.units = 'm/a'

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
  #ncmask.Standard_name = ''
  nccrate.Standard_name = 'calving_rate'

  ncafl.Standard_name  = 'floating_ice_shelf_area'
  ncagr.Standard_name  = 'grounded_ice_sheet_area'
  nclim.Standard_name  =  'land_ice_mass'
  nclimw.Standard_name = 'land_ice_mass_not_displacing_sea_water'
  nctcf.Standard_name  =  'tendency_of_land_ice_mass_due_to_calving'
  nctgf.Standard_name  =  'tendency_of_grounded_ice_mass'

  ncatnw.Standard_name  = 'total_ice_area_NorthWest'
  ncatne.Standard_name  = 'total_ice_area_NorthEast'
  ncatsw.Standard_name  = 'total_ice_area_SouthWest'
  ncatse.Standard_name  = 'total_ice_area_SouthEast'

  ncsa.Standard_name   = 'distance_along_profile_A'
  ncthka.Standard_name = 'land_ice_thickness_along_profile_A'
  ncxvma.Standard_name = 'land_ice_vertical_mean_x_velocity_along_profile_A'
  ncyvma.Standard_name = 'land_ice_vertical_mean_y_velocity_along_profile_A'

  ncxcfa.Standard_name   = 'x_calving_front_on_profile_A'
  ncycfa.Standard_name   = 'y_calving_front_on_profile_A'
  ncxvmcfa.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_profile_A'
  ncyvmcfa.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_profile_A'
  ncthkcfa.Standard_name = 'land_ice_thickness_at_calving_front_on_profile_A'

  ncsb.Standard_name   = 'distance_along_profile_B'
  ncthkb.Standard_name = 'land_ice_thickness_along_profile_B'
  ncxvmb.Standard_name = 'land_ice_vertical_mean_x_velocity_along_profile_B'
  ncyvmb.Standard_name = 'land_ice_vertical_mean_y_velocity_along_profile_B'

  ncxcfb.Standard_name   = 'x_calving_front_on_profile_B'
  ncycfb.Standard_name   = 'y_calving_front_on_profile_B'
  ncxvmcfb.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_profile_B'
  ncyvmcfb.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_profile_B'
  ncthkcfb.Standard_name = 'land_ice_thickness_at_calving_front_on_profile_B'

  ncsc.Standard_name   = 'distance_along_profile_C'
  ncthkc.Standard_name = 'land_ice_thickness_along_profile_C'
  ncxvmc.Standard_name = 'land_ice_vertical_mean_x_velocity_along_profile_C'
  ncyvmc.Standard_name = 'land_ice_vertical_mean_y_velocity_along_profile_C'

  ncxcfc.Standard_name   = 'x_calving_front_on_profile_C'
  ncycfc.Standard_name   = 'y_calving_front_on_profile_C'
  ncxvmcfc.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_profile_C'
  ncyvmcfc.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_profile_C'
  ncthkcfc.Standard_name = 'land_ice_thickness_at_calving_front_on_profile_C'

  ncsd.Standard_name   = 'distance_along_profile_D'
  ncthkd.Standard_name = 'land_ice_thickness_along_profile_D'
  ncxvmd.Standard_name = 'land_ice_vertical_mean_x_velocity_along_profile_D'
  ncyvmd.Standard_name = 'land_ice_vertical_mean_y_velocity_along_profile_D'

  ncxcfd.Standard_name   = 'x_calving_front_on_profile_D'
  ncycfd.Standard_name   = 'y_calving_front_on_profile_D'
  ncxvmcfd.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_profile_D'
  ncyvmcfd.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_profile_D'
  ncthkcfd.Standard_name = 'land_ice_thickness_at_calving_front_on_profile_D'

  ncse.Standard_name   = 'distance_along_profile_E'
  ncthke.Standard_name = 'land_ice_thickness_along_profile_E'
  ncxvme.Standard_name = 'land_ice_vertical_mean_x_velocity_along_profile_E'
  ncyvme.Standard_name = 'land_ice_vertical_mean_y_velocity_along_profile_E'

  ncxcfe.Standard_name   = 'x_calving_front_on_profile_E'
  ncycfe.Standard_name   = 'y_calving_front_on_profile_e'
  ncxvmcfe.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_profile_E'
  ncyvmcfe.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_profile_E'
  ncthkcfe.Standard_name = 'land_ice_thickness_at_calving_front_on_profile_E'

  ncsf.Standard_name   = 'distance_along_profile_F'
  ncthkf.Standard_name = 'land_ice_thickness_along_profile_F'
  ncxvmf.Standard_name = 'land_ice_vertical_mean_x_velocity_along_profile_F'
  ncyvmf.Standard_name = 'land_ice_vertical_mean_y_velocity_along_profile_F'

  ncxcff.Standard_name   = 'x_calving_front_on_profile_F'
  ncycff.Standard_name   = 'y_calving_front_on_profile_F'
  ncxvmcff.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_profile_F'
  ncyvmcff.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_profile_F'
  ncthkcff.Standard_name = 'land_ice_thickness_at_calving_front_on_profile_F'

  ncsg.Standard_name   = 'distance_along_profile_G'
  ncthkg.Standard_name = 'land_ice_thickness_along_profile_G'
  ncxvmg.Standard_name = 'land_ice_vertical_mean_x_velocity_along_profile_G'
  ncyvmg.Standard_name = 'land_ice_vertical_mean_y_velocity_along_profile_G'

  ncxcfg.Standard_name   = 'x_calving_front_on_profile_G'
  ncycfg.Standard_name   = 'y_calving_front_on_profile_G'
  ncxvmcfg.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_profile_G'
  ncyvmcfg.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_profile_G'
  ncthkcfg.Standard_name = 'land_ice_thickness_at_calving_front_on_profile_G'

  ncsh.Standard_name   = 'distance_along_profile_H'
  ncthkh.Standard_name = 'land_ice_thickness_along_profile_H'
  ncxvmh.Standard_name = 'land_ice_vertical_mean_x_velocity_along_profile_H'
  ncyvmh.Standard_name = 'land_ice_vertical_mean_y_velocity_along_profile_H'

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

