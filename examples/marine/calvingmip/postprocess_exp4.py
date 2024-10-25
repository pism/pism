#!/usr/bin/env python3

# Creates output file from PISM EXP4 result to upload for CalvinMIP, 
# as instructed from https://github.com/JRowanJordan/CalvingMIP/wiki/Experiment-4

import numpy as np
import netCDF4 as nc
import datetime
import postprocess_helper as ph

secperyear=365*24*3600

resolution=5.0
#dkm=10.0 #km steps
dkm=5.0 #km steps
vcr=1e-8 # terminal velocity threshold


pismpath     = "/p/tmp/albrecht/pism23/calvmip/thule/exp4-05km-dir-retreat/"
pism_outfile = pismpath + "results/extra_exp4.nc"
pism_tsfile  = pismpath + "results/ts_exp4.nc"

pism_infile = pismpath + "input/thule_input_5km.nc"

exp_outfile = "CalvingMIP_EXP4_PISM_PIK.nc"

#################################################################################################

print('Load PISM data...')
#tex=-1
#crop the outer band of 200km thickness
cr1=40
cr2=361

cr1=0
cr2=401

datnc = nc.Dataset(pism_outfile,"r")
exp_x = datnc.variables["x"][cr1:cr2]
exp_y = datnc.variables["y"][cr1:cr2]
exp_t = datnc.variables["time"][:]/secperyear-1.1e5
try:
    exp_xvm = datnc.variables["xvelmean"][:,cr1:cr2,cr1:cr2]*secperyear
    exp_yvm = datnc.variables["yvelmean"][:,cr1:cr2,cr1:cr2]*secperyear
    exp_thk = datnc.variables["lithk"][:,cr1:cr2,cr1:cr2]
except:
    exp_xvm = datnc.variables["u_ssa"][:,cr1:cr2,cr1:cr2]*secperyear
    exp_yvm = datnc.variables["v_ssa"][:,cr1:cr2,cr1:cr2]*secperyear
    exp_thk = datnc.variables["thk"][:,cr1:cr2,cr1:cr2]
exp_mask = datnc.variables["mask"][:,cr1:cr2,cr1:cr2]-1
exp_crate = datnc.variables["calvingmip_calving_rate"][:,cr1:cr2,cr1:cr2]*secperyear
exp_topg = datnc.variables["topg"][:,cr1:cr2,cr1:cr2]
exp_subh = datnc.variables["ice_area_specific_volume"][:,cr1:cr2,cr1:cr2]
datnc.close()


exp_thks=exp_thk+exp_subh # sum of thk and ice_area_specific_volume
exp_thkm=np.ma.masked_array(exp_thk,mask=exp_thk<=0.0) # exclude empty cells from mean
exp_subhm=np.ma.masked_array(exp_subh,mask=exp_subh<=0.0) # exclude empty subgrid cells

exp_xvmc=np.ma.masked_array(exp_xvm,mask=np.abs(exp_xvm)<=vcr) # remove small values
exp_yvmc=np.ma.masked_array(exp_yvm,mask=np.abs(exp_yvm)<=vcr)


print(np.shape(exp_mask))
#print(exp_t)
Mt,Mx,My = np.shape(exp_mask)

exp_thk[exp_thk == 0.0] = np.nan
exp_xvm[exp_xvm == 0.0] = np.nan
exp_yvm[exp_yvm == 0.0] = np.nan

# avoid empty first snapshot
exp_xvm[0]=exp_xvm[1]
exp_yvm[0]=exp_yvm[1]

datnc = nc.Dataset(pism_tsfile,"r")
exp_ts_t      = datnc.variables["time"][:]/secperyear-1.1e5
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
print(Mx,My,Mp)

MC1=int(390.0/resolution)
MC2=int(590.0/resolution)
MC3=int(450.0/resolution)
MH1=int(150.0/resolution)
MH2=int(740.0/resolution)

#double length
MC2=int(790.0/resolution)
MC3=int(900.0/resolution)
MH2=int(900.0/resolution)

MC2=int(690.0/resolution)
MC3=int(675.0/resolution)
MH2=int(820.0/resolution)

points_CA = [[Mp-MC1,Mp],[Mp-MC2,Mp+MC3]]
points_CB = [[Mp+MC1,Mp],[Mp+MC2,Mp+MC3]]
points_CC = [[Mp-MC1,Mp],[Mp-MC2,Mp-MC3]]
points_CD = [[Mp+MC1,Mp],[Mp+MC2,Mp-MC3]]

points_HA = [[Mp-MH1,Mp],[Mp-MH1,Mp+MH2]]
points_HB = [[Mp+MH1,Mp],[Mp+MH1,Mp+MH2]]
points_HC = [[Mp-MH1,Mp],[Mp-MH1,Mp-MH2]]
points_HD = [[Mp+MH1,Mp],[Mp+MH1,Mp-MH2]]

#switched x and y
#points_CA = [[Mp-MC2,Mp+MC3],[Mp-MC1,Mp]]
#points_CB = [[Mp+MC2,Mp+MC3],[Mp+MC1,Mp]]
#points_CC = [[Mp-MC2,Mp-MC3],[Mp-MC1,Mp]]
#points_CD = [[Mp+MC2,Mp-MC3],[Mp+MC1,Mp]]

#points_HA = [[Mp-MH1,Mp+MH2],[Mp-MH1,Mp]]
#points_HB = [[Mp+MH1,Mp+MH2],[Mp+MH1,Mp]]
#points_HC = [[Mp-MH1,Mp-MH2],[Mp-MH1,Mp]]
#points_HD = [[Mp+MH1,Mp-MH2],[Mp+MH1,Mp]]




transects=[points_CA,points_CB,points_CC,points_CD,points_HA,points_HB,points_HC,points_HD]
point_names=['Caprona A','Caprona B','Caprona C','Caprona D','Halbrane A','Halbrane B','Halbrane C','Halbrane D']

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
  #print(' ')
  for l,po in enumerate(trans):
    profile=[]
    dx=po[1][1]-po[0][1]
    dy=po[1][0]-po[0][0]
    for k,p in enumerate(po):
        
        i=int(np.floor(p[1]))
        j=int(np.floor(p[0]))
        dj=p[0]-np.floor(p[0])
        di=p[1]-np.floor(p[1])

        if ti==0:
            xav[l,k]=ph.interpolate_along_transect(exp_ym,i,j,di,dj)
            yav[l,k]=ph.interpolate_along_transect(exp_xm,i,j,di,dj)
             
            sav[l,k]=np.sqrt((xav[l,k]-xav[l,0])**2+(yav[l,k]-yav[l,0])**2)
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


        if (exp_mask[ti,int(np.sign(dx)+i),int(np.sign(dy)+j)]-exp_mask[ti,i,j]==1.0):

            Atot    = 9.0 # area of neighbor cells around i,j
            thkmean = np.nanmean(exp_thkm[ti,i-1:i+2,j-1:j+2]) # mean over cells with thk>0
            Acov    = np.sum((exp_thks[ti,i-1:i+2,j-1:j+2])/thkmean) # area covered by cells with thk>0 or icea_area_specific_volume>0
            ds      = 3.0*(Acov/Atot-0.5) # average distance of i,j from ice front, assuming half of 9 cells would be filled with thk>0 
            ycf[ti,l]=exp_ym[i,j]+ds*dy*(exp_ym[i,j+1]-exp_ym[i,j])
            xcf[ti,l]=exp_xm[i,j]+ds*dx*(exp_xm[i+1,j]-exp_xm[i,j])

            Hcf[ti,l]=thkmean
            try:
                vxcf[ti,l]=np.nanmean(exp_xvmc[ti,i-1:i+2,j-1:j+2])
                vycf[ti,l]=np.nanmean(exp_yvmc[ti,i-1:i+2,j-1:j+2])
            except:
                vxcf[ti,l]=vxcf[ti-1,l]
                vycf[ti,l]=vycf[ti-1,l]

            #print(ti,l,i,j,xcf[ti,l],ycf[ti,l],Hcf[ti,l])


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


######################################################################################################
#cut off outer region

cr1=40
cr2=361


########################################################################################################
lent=Mt

if True:
  print("Write data to netCDF file "+exp_outfile)
  wrtfile = nc.Dataset(exp_outfile, 'w', format='NETCDF4_CLASSIC')

  wrtfile.createDimension('X', size=len(exp_x[cr1:cr2]))
  wrtfile.createDimension('Y', size=len(exp_y[cr1:cr2]))

  wrtfile.createDimension('Time1', size=len(exp_t)) # every year
  wrtfile.createDimension('Time100', size=len(exp_t[::100])) # every 100 years


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
    
  wrtfile.createDimension('Caprona A', size=np.shape(trans[0])[0])
  ncxa     = wrtfile.createVariable('Caprona A', 'f8', ('Caprona A',))
  ncsa    = wrtfile.createVariable('sCapA', 'f8', ('Caprona A'))
  ncthka  = wrtfile.createVariable('lithkCapA', 'f8', ('Time1','Caprona A'),fill_value=np.nan)
  ncxvma  = wrtfile.createVariable('xvelmeanCapA', 'f8', ('Time1','Caprona A'),fill_value=np.nan)
  ncyvma  = wrtfile.createVariable('yvelmeanCapA', 'f8', ('Time1','Caprona A'),fill_value=np.nan)
  ncmaska  = wrtfile.createVariable('maskCapA', 'f8', ('Time1','Caprona A'))
    
  wrtfile.createDimension('Caprona B', size=np.shape(trans[1])[0])
  ncxb     = wrtfile.createVariable('Caprona B', 'f8', ('Caprona B',))
  ncsb    = wrtfile.createVariable('sCapB', 'f8', ('Caprona B'))
  ncthkb  = wrtfile.createVariable('lithkCapB', 'f8', ('Time1','Caprona B'),fill_value=np.nan)
  ncxvmb  = wrtfile.createVariable('xvelmeanCapB', 'f8', ('Time1','Caprona B'),fill_value=np.nan)
  ncyvmb  = wrtfile.createVariable('yvelmeanCapB', 'f8', ('Time1','Caprona B'),fill_value=np.nan)
  ncmaskb  = wrtfile.createVariable('maskCapB', 'f8', ('Time1','Caprona B'))

  wrtfile.createDimension('Caprona C', size=np.shape(trans[2])[0])
  ncxc     = wrtfile.createVariable('Caprona C', 'f8', ('Caprona C',))
  ncsc    = wrtfile.createVariable('sCapC', 'f8', ('Caprona C'))
  ncthkc  = wrtfile.createVariable('lithkCapC', 'f8', ('Time1','Caprona C'),fill_value=np.nan)
  ncxvmc  = wrtfile.createVariable('xvelmeanCapC', 'f8', ('Time1','Caprona C'),fill_value=np.nan)
  ncyvmc  = wrtfile.createVariable('yvelmeanCapC', 'f8', ('Time1','Caprona C'),fill_value=np.nan)
  ncmaskc  = wrtfile.createVariable('maskCapC', 'f8', ('Time1','Caprona C'))
    
  wrtfile.createDimension('Caprona D', size=np.shape(trans[3])[0])
  ncxd     = wrtfile.createVariable('Caprona D', 'f8', ('Caprona D',))
  ncsd    = wrtfile.createVariable('sCapD', 'f8', ('Caprona D'))
  ncthkd  = wrtfile.createVariable('lithkCapD', 'f8', ('Time1','Caprona D'),fill_value=np.nan)
  ncxvmd  = wrtfile.createVariable('xvelmeanCapD', 'f8', ('Time1','Caprona D'),fill_value=np.nan)
  ncyvmd  = wrtfile.createVariable('yvelmeanCapD', 'f8', ('Time1','Caprona D'),fill_value=np.nan)
  ncmaskd  = wrtfile.createVariable('maskCapD', 'f8', ('Time1','Caprona D'))

  wrtfile.createDimension('Halbrane A', size=np.shape(trans[4])[0])
  ncxe     = wrtfile.createVariable('Halbrane A', 'f8', ('Halbrane A',))
  ncse    = wrtfile.createVariable('sHalA', 'f8', ('Halbrane A'))
  ncthke  = wrtfile.createVariable('lithkHalA', 'f8', ('Time1','Halbrane A'),fill_value=np.nan)
  ncxvme  = wrtfile.createVariable('xvelmeanHalA', 'f8', ('Time1','Halbrane A'),fill_value=np.nan)
  ncyvme  = wrtfile.createVariable('yvelmeanHalA', 'f8', ('Time1','Halbrane A'),fill_value=np.nan)
  ncmaske  = wrtfile.createVariable('maskHalA', 'f8', ('Time1','Halbrane A'))
    
  wrtfile.createDimension('Halbrane B', size=np.shape(trans[5])[0])
  ncxf     = wrtfile.createVariable('Halbrane B', 'f8', ('Halbrane B',))
  ncsf    = wrtfile.createVariable('sHalB', 'f8', ('Halbrane B'))
  ncthkf  = wrtfile.createVariable('lithkHalB', 'f8', ('Time1','Halbrane B'),fill_value=np.nan)
  ncxvmf  = wrtfile.createVariable('xvelmeanHalB', 'f8', ('Time1','Halbrane B'),fill_value=np.nan)
  ncyvmf  = wrtfile.createVariable('yvelmeanHalB', 'f8', ('Time1','Halbrane B'),fill_value=np.nan)
  ncmaskf  = wrtfile.createVariable('maskHalB', 'f8', ('Time1','Halbrane B'))

  wrtfile.createDimension('Halbrane C', size=np.shape(trans[6])[0])
  ncxg     = wrtfile.createVariable('Halbrane C', 'f8', ('Halbrane C',))
  ncsg    = wrtfile.createVariable('sHalC', 'f8', ('Halbrane C'))
  ncthkg  = wrtfile.createVariable('lithkHalC', 'f8', ('Time1','Halbrane C'),fill_value=np.nan)
  ncxvmg  = wrtfile.createVariable('xvelmeanHalC', 'f8', ('Time1','Halbrane C'),fill_value=np.nan)
  ncyvmg  = wrtfile.createVariable('yvelmeanHalC', 'f8', ('Time1','Halbrane C'),fill_value=np.nan)
  ncmaskg  = wrtfile.createVariable('maskHalC', 'f8', ('Time1','Halbrane C'))
    
  wrtfile.createDimension('Halbrane D', size=np.shape(trans[7])[0])
  ncxh     = wrtfile.createVariable('Halbrane D', 'f8', ('Halbrane D',))
  ncsh    = wrtfile.createVariable('sHalD', 'f8', ('Halbrane D'))
  ncthkh  = wrtfile.createVariable('lithkHalD', 'f8', ('Time1','Halbrane D'),fill_value=np.nan)
  ncxvmh  = wrtfile.createVariable('xvelmeanHalD', 'f8', ('Time1','Halbrane D'),fill_value=np.nan)
  ncyvmh  = wrtfile.createVariable('yvelmeanHalD', 'f8', ('Time1','Halbrane D'),fill_value=np.nan)
  ncmaskh  = wrtfile.createVariable('maskHalD', 'f8', ('Time1','Halbrane D'))


  ncatnw  = wrtfile.createVariable('iareatotalNW', 'f8', ('Time1'))
  ncatne  = wrtfile.createVariable('iareatotalNE', 'f8', ('Time1'))
  ncatsw  = wrtfile.createVariable('iareatotalSW', 'f8', ('Time1'))
  ncatse  = wrtfile.createVariable('iareatotalSE', 'f8', ('Time1'))

  ncxcfa  = wrtfile.createVariable('xcfCapA', 'f8', ('Time1'),fill_value=np.nan)
  ncycfa  = wrtfile.createVariable('ycfCapA', 'f8', ('Time1'),fill_value=np.nan)
  ncxvmcfa  = wrtfile.createVariable('xvelmeancfCapA', 'f8', ('Time1'))
  ncyvmcfa  = wrtfile.createVariable('yvelmeancfCapA', 'f8', ('Time1'))
  ncthkcfa  = wrtfile.createVariable('lithkcfCapA', 'f8', ('Time1'))

  ncxcfb  = wrtfile.createVariable('xcfCapB', 'f8', ('Time1'),fill_value=np.nan)
  ncycfb  = wrtfile.createVariable('ycfCapB', 'f8', ('Time1'),fill_value=np.nan)
  ncxvmcfb  = wrtfile.createVariable('xvelmeancfCapB', 'f8', ('Time1'))
  ncyvmcfb  = wrtfile.createVariable('yvelmeancfCapB', 'f8', ('Time1'))
  ncthkcfb  = wrtfile.createVariable('lithkcfCapB', 'f8', ('Time1'))

  ncxcfc  = wrtfile.createVariable('xcfCapC', 'f8', ('Time1'),fill_value=np.nan)
  ncycfc  = wrtfile.createVariable('ycfCapC', 'f8', ('Time1'),fill_value=np.nan)
  ncxvmcfc  = wrtfile.createVariable('xvelmeancfCapC', 'f8', ('Time1'))
  ncyvmcfc  = wrtfile.createVariable('yvelmeancfCapC', 'f8', ('Time1'))
  ncthkcfc  = wrtfile.createVariable('lithkcfCapC', 'f8', ('Time1'))

  ncxcfd  = wrtfile.createVariable('xcfCapD', 'f8', ('Time1'),fill_value=np.nan)
  ncycfd  = wrtfile.createVariable('ycfCapD', 'f8', ('Time1'),fill_value=np.nan)
  ncxvmcfd  = wrtfile.createVariable('xvelmeancfCapD', 'f8', ('Time1'))
  ncyvmcfd  = wrtfile.createVariable('yvelmeancfCapD', 'f8', ('Time1'))
  ncthkcfd  = wrtfile.createVariable('lithkcfCapD', 'f8', ('Time1'))

  ncxcfe  = wrtfile.createVariable('xcfHalA', 'f8', ('Time1'), fill_value=np.nan)
  ncycfe  = wrtfile.createVariable('ycfHalA', 'f8', ('Time1'),fill_value=np.nan)
  ncxvmcfe  = wrtfile.createVariable('xvelmeancfHalA', 'f8', ('Time1'))
  ncyvmcfe  = wrtfile.createVariable('yvelmeancfHalA', 'f8', ('Time1'))
  ncthkcfe  = wrtfile.createVariable('lithkcfHalA', 'f8', ('Time1'))

  ncxcff  = wrtfile.createVariable('xcfHalB', 'f8', ('Time1'),fill_value=np.nan)
  ncycff  = wrtfile.createVariable('ycfHalB', 'f8', ('Time1'),fill_value=np.nan)
  ncxvmcff  = wrtfile.createVariable('xvelmeancfHalB', 'f8', ('Time1'))
  ncyvmcff  = wrtfile.createVariable('yvelmeancfHalB', 'f8', ('Time1'))
  ncthkcff  = wrtfile.createVariable('lithkcfHalB', 'f8', ('Time1'))

  ncxcfg  = wrtfile.createVariable('xcfHalC', 'f8', ('Time1'),fill_value=np.nan)
  ncycfg  = wrtfile.createVariable('ycfHalC', 'f8', ('Time1'),fill_value=np.nan)
  ncxvmcfg  = wrtfile.createVariable('xvelmeancfHalC', 'f8', ('Time1'))
  ncyvmcfg  = wrtfile.createVariable('yvelmeancfHalC', 'f8', ('Time1'))
  ncthkcfg  = wrtfile.createVariable('lithkcfHalC', 'f8', ('Time1'))

  ncxcfh  = wrtfile.createVariable('xcfHalD', 'f8', ('Time1'),fill_value=np.nan)
  ncycfh  = wrtfile.createVariable('ycfHalD', 'f8', ('Time1'),fill_value=np.nan)
  ncxvmcfh  = wrtfile.createVariable('xvelmeancfHalD', 'f8', ('Time1'))
  ncyvmcfh  = wrtfile.createVariable('yvelmeancfHalD', 'f8', ('Time1'))
  ncthkcfh  = wrtfile.createVariable('lithkcfHalD', 'f8', ('Time1'))
    
  #####################################################################################
    
  nct[:] = exp_t[:]
  nct2[:]= exp_t[::100]


  ncx[:] = exp_x[cr1:cr2]
  ncy[:] = exp_y[cr1:cr2]
    
  ncxvm[:]  = exp_xvm[::100,cr1:cr2,cr1:cr2]
  ncyvm[:]  = exp_yvm[::100,cr1:cr2,cr1:cr2]
  ncthk[:]  = exp_thk[::100,cr1:cr2,cr1:cr2]
  ncmask[:] = exp_mask[::100,cr1:cr2,cr1:cr2]
  nccrate[:]= exp_crate[::100,cr1:cr2,cr1:cr2]
  nctopg[:] = exp_topg[::100,cr1:cr2,cr1:cr2]


  #ncx[:] = exp_x_cr[:]
  #ncy[:] = exp_y_cr[:]

  #ncxvm[:]  = exp_xvm_cr[::100]
  #ncyvm[:]  = exp_yvm_cr[::100]
  #ncthk[:]  = exp_thk_cr[::100]
  #ncmask[:] = exp_mask_cr[::100]
  #nccrate[:]= exp_crate_cr[::100]
  #nctopg[:] = exp_topg_cr[::100]

    
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
  

  
  ncsa.Standard_name   = 'distance_along_Caprona_A'
  ncthka.Standard_name = 'land_ice_thickness_along_Caprona_A'
  ncxvma.Standard_name = 'land_ice_vertical_mean_x_velocity_along_Caprona_A'
  ncyvma.Standard_name = 'land_ice_vertical_mean_y_velocity_along_Caprona_A'
    
  ncsb.Standard_name   = 'distance_along_Caprona_B'
  ncthkb.Standard_name = 'land_ice_thickness_along_Caprona_B'
  ncxvmb.Standard_name = 'land_ice_vertical_mean_x_velocity_along_Caprona_B'
  ncyvmb.Standard_name = 'land_ice_vertical_mean_y_velocity_along_Caprona_B'    
    
  ncsc.Standard_name   = 'distance_along_Caprona_C'
  ncthkc.Standard_name = 'land_ice_thickness_along_Caprona_C'
  ncxvmc.Standard_name = 'land_ice_vertical_mean_x_velocity_along_Caprona_C'
  ncyvmc.Standard_name = 'land_ice_vertical_mean_y_velocity_along_Caprona_C'
    
  ncsd.Standard_name   = 'distance_along_Caprona_D'
  ncthkd.Standard_name = 'land_ice_thickness_along_Caprona_D'
  ncxvmd.Standard_name = 'land_ice_vertical_mean_x_velocity_along_Caprona_D'
  ncyvmd.Standard_name = 'land_ice_vertical_mean_y_velocity_along_Caprona_D'
    
  ncse.Standard_name   = 'distance_along_Halbrane_A'
  ncthke.Standard_name = 'land_ice_thickness_along_Halbrane_A'
  ncxvme.Standard_name = 'land_ice_vertical_mean_x_velocity_along_Halbrane_A'
  ncyvme.Standard_name = 'land_ice_vertical_mean_y_velocity_along_Halbrane_A'
    
  ncsf.Standard_name   = 'distance_along_Halbrane_B'
  ncthkf.Standard_name = 'land_ice_thickness_along_Halbrane_B'
  ncxvmf.Standard_name = 'land_ice_vertical_mean_x_velocity_along_Halbrane_B'
  ncyvmf.Standard_name = 'land_ice_vertical_mean_y_velocity_along_Halbrane_B'
    
  ncsg.Standard_name   = 'distance_along_Halbrane_C'
  ncthkg.Standard_name = 'land_ice_thickness_along_Halbrane_C'
  ncxvmg.Standard_name = 'land_ice_vertical_mean_x_velocity_along_Halbrane_C'
  ncyvmg.Standard_name = 'land_ice_vertical_mean_y_velocity_along_Halbrane_C'
    
  ncsh.Standard_name   = 'distance_along_Halbrane_D'
  ncthkh.Standard_name = 'land_ice_thickness_along_Halbrane_D'
  ncxvmh.Standard_name = 'land_ice_vertical_mean_x_velocity_along_Halbrane_D'
  ncyvmh.Standard_name = 'land_ice_vertical_mean_y_velocity_along_Halbrane_D'
    
  ncatnw.Standard_name  = 'total_ice_area_NorthWest'
  ncatne.Standard_name  = 'total_ice_area_NorthEast'
  ncatsw.Standard_name  = 'total_ice_area_SouthWest'
  ncatse.Standard_name  = 'total_ice_area_SouthEast'

  ncxcfa.Standard_name   = 'x_calving_front_on_Caprona_A'
  ncycfa.Standard_name   = 'y_calving_front_on_Caprona_A'
  ncxvmcfa.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_Caprona_A'
  ncyvmcfa.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_Caprona_A'
  ncthkcfa.Standard_name = 'land_ice_thickness_at_calving_front_on_Caprona_A'

  ncxcfb.Standard_name   = 'x_calving_front_on_Caprona_B'
  ncycfb.Standard_name   = 'y_calving_front_on_Caprona_B'
  ncxvmcfb.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_Caprona_B'
  ncyvmcfb.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_Caprona_B'
  ncthkcfb.Standard_name = 'land_ice_thickness_at_calving_front_on_Caprona_B'

  ncxcfc.Standard_name   = 'x_calving_front_on_Caprona_C'
  ncycfc.Standard_name   = 'y_calving_front_on_Caprona_C'
  ncxvmcfc.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_Caprona_C'
  ncyvmcfc.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_Caprona_C'
  ncthkcfc.Standard_name = 'land_ice_thickness_at_calving_front_on_Caprona_C'

  ncxcfd.Standard_name   = 'x_calving_front_on_Caprona_D'
  ncycfd.Standard_name   = 'y_calving_front_on_Caprona_D'
  ncxvmcfd.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_Caprona_D'
  ncyvmcfd.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_Caprona_D'
  ncthkcfd.Standard_name = 'land_ice_thickness_at_calving_front_on_Caprona_D'

  ncxcfe.Standard_name   = 'x_calving_front_on_Halbrane_A'
  ncycfe.Standard_name   = 'y_calving_front_on_Halbrane_A'
  ncxvmcfe.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_Halbrane_A'
  ncyvmcfe.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_Halbrane_A'
  ncthkcfe.Standard_name = 'land_ice_thickness_at_calving_front_on_Halbrane_A'

  ncxcff.Standard_name   = 'x_calving_front_on_Halbrane_B'
  ncycff.Standard_name   = 'y_calving_front_on_Halbrane_B'
  ncxvmcff.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_Halbrane_B'
  ncyvmcff.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_Halbrane_B'
  ncthkcff.Standard_name = 'land_ice_thickness_at_calving_front_on_Halbrane_B'

  ncxcfg.Standard_name   = 'x_calving_front_on_Halbrane_C'
  ncycfg.Standard_name   = 'y_calving_front_on_Halbrane_C'
  ncxvmcfg.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_Halbrane_C'
  ncyvmcfg.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_Halbrane_C'
  ncthkcfg.Standard_name = 'land_ice_thickness_at_calving_front_on_Halbrane_C'

  ncxcfh.Standard_name   = 'x_calving_front_on_Halbrane_D'
  ncycfh.Standard_name   = 'y_calving_front_on_Halbrane_D'
  ncxvmcfh.Standard_name = 'land_ice_vertical_mean_x_velocity_at_calving_front_on_Halbrane_D'
  ncyvmcfh.Standard_name = 'land_ice_vertical_mean_y_velocity_at_calving_front_on_Halbrane_D'
  ncthkcfh.Standard_name = 'land_ice_thickness_at_calving_front_on_Halbrane_D'



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

