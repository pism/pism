#!/usr/bin/env python
from pylab import *
import os
import sys
from getopt import getopt, GetoptError
from netCDF import Dataset as NC

# description of where the transcripts are; EXAMPLES!!
# (0,       1,        2,      3,       4,      5,      6)
# (FILENAME,STARTLINE,ENDLINE,EXPERNUM,GRIDRES,[0|1|2],[P1contSTAGE])
desc = [
('../pst_P3etc.out',34216,43806,1,15,0),  # 15km (coarse) grids
('../pst_P3etc.out',43808,53326,2,15,0),
('../pst_P3etc.out',53328,67690,3,15,0),
('../pst_P3etc.out',67692,77660,4,15,0),
('../pst.out',186263,207682,1,10,0),  # 10km grids
('../pst.out',207684,229045,2,10,0),
('../pst_P3etc.out',2,34214,3,10,0),
('../pst.out',229047,251562,4,10,0),
('../pst_P3etc.out',260878,282524,1,10,1),  # 10km grids w vert refine
('../pst_P3etc.out',282526,304142,2,10,1),
('../pst_P3etc.out',304144,336238,3,10,1),
('../pst_P3etc.out',336240,358826,4,10,1),
('../pst_P3etc.out',77662,115135,1,7.5,0),  # 7.5km (fine) grids
('../pst_P3etc.out',115137,154098,2,7.5,0),
('../pst_P3etc.out',154100,220091,3,7.5,0),
('../pst_P3etc.out',220093,260876,4,7.5,0),
('../mid_6july/pstfinish.pbs.o510189',4,87292,1,5,0),  # 5km (finest) grids
('../mid_6july/pstfinish.pbs.o510189',87296,174820,2,5,0),
('../mid_6july/pstfinish.pbs.o510189',174823,315350,3,5,0),
('../mid_8july/pstfinish2.pbs.o511587',4,86959,4,5,0),
('../mid_7july/pstP1cont.pbs.o510058',4,66949,1,10,2,20),  # P1cont
('../mid_7july/pstP1cont.pbs.o510058',66953,159004,1,10,2,40),
('../mid_7july/pstP1cont.pbs.o510058',159008,252293,1,10,2,60),
('../mid_7july/pstP1cont2.pbs.o511218',4,94013,1,10,2,80),
('../mid_7july/pstP1cont2.pbs.o511218',94017,188112,1,10,2,100),
#('../pst_P1cont.out',2,67385,1,10,2,20),  # P1cont *with -tempskip 10*
#('../pst_P1cont.out',67387,161674,1,10,2,40),
#('../pst_P1cont.out',161676,258755,1,10,2,60),
#('../pst_P1cont.out',258757,357322,1,10,2,80),
#('../pst_P1cont.out',357324,457109,1,10,2,100),
]


colors = ["red", "green", "yellow", "blue", "black", "cyan", "magenta"]
def getcolor(gridspacing):
  if (gridspacing == 5):
    return colors[0]
  elif (gridspacing == 7.5):
    return colors[1]
  elif (gridspacing == 10):
    return colors[2]
  elif (gridspacing == 15):
    return colors[3]
  else:
    print "INVALID HORIZONTAL GRID SPACING; ENDING"
    sys.exit(2)

volstyles = ['-','-.',':','--']
def getvolstyle(gridspacing):
  if (gridspacing == 5):
    return volstyles[0]
  elif (gridspacing == 7.5):
    return volstyles[1]
  elif (gridspacing == 10):
    return volstyles[2]
  elif (gridspacing == 15):
    return volstyles[3]
  else:
    print "INVALID HORIZONTAL GRID SPACING; ENDING"
    sys.exit(2)

P1widths = [70, 30, 100, 50]
P1styles = ['-',':','-.','--']
P2angles = [0,10,45]
P2styles = ['-',':','--']
#P3drops = [0, 500, 1000, 2000]
P3slopepercents = [0.000, 0.077, 0.154, 0.308]
P4phis = [2, 3, 8, 10]

P1contlines = []

numfigs=12

doExtract = False
try:
  opts, args = getopt(sys.argv[1:], "e:",["extract="])
  for opt, arg in opts:
    if opt in ("-e", "--extract"):
      doExtract = True
except GetoptError:
  print "Incorrect command line arguments. Exiting..."
  sys.exit(-1)


print "creating empty figures 1,..,%d" % numfigs
for j in range(1,numfigs+1):
  figure(j)

for k in desc:
  pre = "P%d_%dkm" % (k[3],k[4])
  if (k[5] == 1):
    pre += "VR"
  elif (k[5] == 2):
    pre += "_%dk" % k[6]
  print ""
  print "POST-PROCESSING STDOUT FROM EXPERIMENT %s:" % pre

  if (doExtract):
    extract = "cat %s | sed -n '%d,%dp' > %s.out" % (k[0],k[1],k[2],pre)
    print "extracting standard out: '%s'" % extract
    try:
      status = os.system(extract)
    except KeyboardInterrupt:  sys.exit(2)
    if status:  sys.exit(status)

    outname = "%s.ser.nc" % pre
    print "using series.py to create NetCDF time series"
    print " ------------------ series.py OUTPUT:"
    try:
      status = os.system("series.py -f %s.out -o %s" % (pre,outname))
    except KeyboardInterrupt:  sys.exit(2)
    if status:  sys.exit(status)
    print " ------------------"

  filename = "%s.ser.nc" % pre
  print "opening %s to make time series figures" % filename
  nc = NC(filename, 'r')
  time = nc.variables["t"][:]

  # figure(1|2) shows grid refinement for ivol time series
  #   for P[1|2]
  if ( (k[3] < 3) & ((k[5] == 0) | (k[5] == 1)) ):
    print "adding ivol to figure(%d)" % (k[3])
    figure(k[3])
    hold(True)
    var = nc.variables["ivol"][:]
    myls = getvolstyle(k[4])
    mylabel = "%.1f km grid" % k[4]
    if (k[5] == 1):
      mylabel += " (fine vert)"
      plot(time, var, linestyle=myls, linewidth=4, color='black', label=mylabel)
    else:
      plot(time, var, linestyle=myls, linewidth=1.5, color='black', label=mylabel)
    hold(False)

  # figure(5) shows P2 down stream speeds from 5km run
  if ( (k[3] == 2) & (k[4] == 5) ):  # want only finest
    print "adding avdwn0,avdwn1,avdwn2 to figure(5)"
    figure(5)
    hold(True)
    for j in [0,1,2]:
      varname = "avdwn%d" % j
      var = nc.variables[varname][:]
      mylabel = r"$%d^\circ$ strip" % P2angles[j]
      semilogy(time, var, linestyle=P2styles[j], linewidth=1.5, 
               color='black', label=mylabel)

  # figure(6) shows P1cont ivol time series
  # figure(7) shows P1cont down stream speed time series
  #     for 30,50,70 km wide streams
  # figure(8) shows P1cont down stream speed time series
  #     for 100km wide stream
  if ( (k[5] == 2) | ((k[3] == 1) & (k[4] == 10) & (k[5] == 0)) ):
    print "adding ivol from P1cont to figure(6)"
    figure(6)
    hold(True)
    var = nc.variables["ivol"][:]
    plot(time, var, linewidth=1.5, color='black')
    hold(False)
    print "adding avdwn013 (70,30,50km wide) from P1cont to figure(7)"
    figure(7)
    hold(True)
    for j in [0,1,3]:
      varname = "avdwn%d" % j
      var = nc.variables[varname][:]
      myl = semilogy(time, var, linestyle=P1styles[j], linewidth=1.5, 
                     color='black')
      P1contlines.append(myl)
    print "adding avdwn2 (100km wide) from P1cont to figure(8)"
    figure(8)
    hold(True)
    varname = "avdwn2"
    var = nc.variables[varname][:]
    myl = semilogy(time, var, linewidth=1.5, color='black')
    P1contlines.append(myl)
    if (k[5] == 2):
      if (k[6] == 100): # if last part, then extract last 5ka
        time_100km_last5ka = time[time>=95000]
        avdwn_100km_last5ka = var[time>=95000]

  # figure(9) shows ivol time series for P[1|2|3|4] on 5km grid
  if ((k[4] == 5) & (k[5] == 0)):
    print "adding ivol to figure(9)"
    figure(9)
    hold(True)
    var = nc.variables["ivol"][:]
    mylabel = "P%d" % k[3]
    plot(time, var, linewidth=1.5, color="black",
         linestyle=P1styles[k[3]-1], label=mylabel)
    hold(False)
    
  # figure(10) shows P1 downstream speeds from 5km run
  if ( (k[3] == 1) & (k[4] == 5) ):  # want only finest
    print "adding avdwn0,avdwn1,avdwn2 to figure(10)"
    figure(10)
    hold(True)
    for j in [0,1,2,3]:
      varname = "avdwn%d" % j
      var = nc.variables[varname][:]
      mylabel = r"$%d$ km" % P1widths[j]
      semilogy(time, var, linestyle=P1styles[j], linewidth=1.5, 
               color='black', label=mylabel)
    
  # figure(11) shows P3 downstream speeds from 5km run
  if ( (k[3] == 3) & (k[4] == 5) ):  # want only finest
    print "adding avdwn0,avdwn1,avdwn2 to figure(11)"
    figure(11)
    hold(True)
    for j in [0,1,2,3]:
      varname = "avdwn%d" % j
      var = nc.variables[varname][:]
      #mylabel = r"$%d$ m" % P3drops[j]
      mylabel = r"$%.3f$ %%" % P3slopepercents[j]
      semilogy(time, var, linestyle=P1styles[j], linewidth=1.5, 
               color='black', label=mylabel)
    
  # figure(12) shows P4 downstream speeds from 5km run
  if ( (k[3] == 4) & (k[4] == 5) ):  # want only finest
    print "adding avdwn0,avdwn1,avdwn2 to figure(12)"
    figure(12)
    hold(True)
    for j in [0,1,2,3]:
      varname = "avdwn%d" % j
      var = nc.variables[varname][:]
      mylabel = r"$%d^\circ$" % P4phis[j]
      semilogy(time, var, linestyle=P1styles[j], linewidth=1.5, 
               color='black', label=mylabel)

  nc.close()
   
  
print ""
print "LABELING AND SAVING FIGURES (P*.png)"

def savepng(fignum,prefix):
  filename = prefix + ".png"
  print "saving figure(%d) as %s ..." % (fignum,filename)
  savefig(filename, dpi=300, facecolor='w', edgecolor='w')
    
# P1 and P2 ivol figures (for grid refinement)
for j in [1,2]:
  figure(j)
  xlabel("t  (a)",fontsize=16)
  ylabel(r'volume  ($10^6$ $\mathrm{km}^3$)',fontsize=16)
  axis([0, 5000, 2.14, 2.22]) 
  xticks(arange(0,6000,1000),fontsize=14)
  yticks(arange(2.14,2.22,0.02),fontsize=14)
  legend(loc='upper right')
  savepng(j,"P%d_vol" % j)

# P2 speed figure
figure(5)
xlabel("t  (a)",fontsize=16)
ylabel("average downstream speed  (m/a)",fontsize=16)
axis([0, 5000, 10, 1000]) 
xticks(arange(0,6000,1000),fontsize=14)
yticks([10,20,50,100,200,500,1000],('10','20','50','100','200','500','1000'),fontsize=14)
rcParams.update({'legend.fontsize': 14})
legend()
savepng(5,"P2_dwnspeeds")
hold(False)

# P1cont volume figure
figure(6)
xlabel("t  (ka)",fontsize=16)
ylabel(r'volume  ($10^6$ $\mathrm{km}^3$)',fontsize=16)
axis([0, 100000, 2.10, 2.22]) 
xticks(arange(0,120000,20000),('0','20','40','60','80','100'),fontsize=14)
yticks(arange(2.10,2.22,0.02),fontsize=14)
savepng(6,"P1cont_vol")
hold(False)

# P1cont speed figure for 70,30,50km wide
figure(7)
xlabel("t  (ka)",fontsize=16)
ylabel("average downstream speed  (m/a)",fontsize=16)
axis([0, 100000, 75, 200]) 
xticks(arange(0,120000,20000),('0','20','40','60','80','100'),fontsize=14)
yticks([75,100,150,200],('75','100','150','200'),fontsize=14)
rcParams.update({'legend.fontsize': 14})
legend(P1contlines[0:3],
       ( "%d km wide strip" % P1widths[0],"%d km wide strip" % P1widths[1],
         "%d km wide strip" % P1widths[3] ))
savepng(7,"P1cont_dwnspeeds")
hold(False)

# P1cont speed figure for 100km wide, with inset
figure(8)
xlabel("t  (ka)",fontsize=16)
ylabel("average downstream speed  (m/a)",fontsize=16)
axis([0, 100000, 75, 200]) 
xticks(arange(0,120000,20000),('0','20','40','60','80','100'),fontsize=14)
yticks([75,100,150,200],('75','100','150','200'),fontsize=14)
hold(False)
axes100km = axes([0.67, 0.71, 0.28, 0.24],axisbg='w')
plot(time_100km_last5ka,avdwn_100km_last5ka, linewidth=1.5, color='black')
setp(axes100km,xticks=arange(95.e3,101.e3,1.e3),
     xticklabels=('95','96','97','98','99','100'))
savepng(8,"P1cont_dwnspeed100km")
    
# ivol figure for all P? at 5km resolution
figure(9)
xlabel("t  (a)",fontsize=16)
ylabel(r'volume  ($10^6$ $\mathrm{km}^3$)',fontsize=16)
axis([0, 5000, 2.14, 2.22]) 
xticks(arange(0,6000,1000),fontsize=14)
yticks(arange(2.14,2.22,0.02),fontsize=14)
legend(loc='upper right')
savepng(9,"Pall_vol")

# P1 speed figure
figure(10)
xlabel("t  (a)",fontsize=16)
ylabel("average downstream speed  (m/a)",fontsize=16)
axis([0, 5000, 10, 1000]) 
xticks(arange(0,6000,1000),fontsize=14)
yticks([10,20,50,100,200,500,1000],('10','20','50','100','200','500','1000'),fontsize=14)
rcParams.update({'legend.fontsize': 14})
legend()
savepng(10,"P1_dwnspeeds")

# P3 speed figure
figure(11)
xlabel("t  (a)",fontsize=16)
ylabel("average downstream speed  (m/a)",fontsize=16)
axis([0, 5000, 10, 1000]) 
xticks(arange(0,6000,1000),fontsize=14)
yticks([10,20,50,100,200,500,1000],('10','20','50','100','200','500','1000'),fontsize=14)
rcParams.update({'legend.fontsize': 14})
legend()
savepng(11,"P3_dwnspeeds")

# P4 speed figure
figure(12)
xlabel("t  (a)",fontsize=16)
ylabel("average downstream speed  (m/a)",fontsize=16)
axis([0, 5000, 10, 1000]) 
xticks(arange(0,6000,1000),fontsize=14)
yticks([10,20,50,100,200,500,1000],('10','20','50','100','200','500','1000'),fontsize=14)
rcParams.update({'legend.fontsize': 14})
legend()
savepng(12,"P4_dwnspeeds")


print ""
print "AUTOCROPPING FIGURES (P*.png)"
# uses one of the ImageMagick tools (http://www.imagemagick.org/)
try:
  status = os.system("mogrify -verbose -trim +repage P*.png")
except KeyboardInterrupt:  sys.exit(2)
if status:  sys.exit(status)


exit()

