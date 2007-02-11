% ROSS_PLOT Plots computed speed of Ross and shows observed values from RIGGS
% data.  Assumes foo.m result of PISM Ross computation has already been run
% (i.e. "$ pisms -ross -o foo -of m" and then do ">> foo"), so variables 
% c,u,v,H,mask,x,y are defined.  Reads RIGGS data from riggs_ELBclean.dat.
% Compare figures 3 from (Bently 1984) and figures 11 and 14 from 
% (Thomas et al 1984) and figures in (MacAyeal et al 1996).

% ELB 2/4/07; 2/11/07

% see 111by147.dat for these ranges
dlat = (-5.42445 - (-12.3325))/110;
gridlatext = linspace(-12.3325 - dlat * 46,-5.42445,147);
gridlon = linspace(-5.26168,3.72207,147);

% load RIGGS data FROM D. MACAYEAL TO ELB ON 19 DEC 2006.
load -ascii riggs_ELBclean.dat
RIGGS=riggs_ELBclean;
clear riggs_ELBclean;

% show velocities as color
cforplot=c;  cforplot(H<20) = -20; cforplot(mask==1) = 1200;
imagesc(gridlon,gridlatext,cforplot'), colorbar
h=get(gcf,'CurrentAxes');  set(h, 'YDir', 'normal')
axis([-5.26168 3.72207 -13 -5.42445])
xlabel('RIGGS grid longitude (deg E)'), ylabel('RIGGS grid latitude (deg N)')

% compute grid lat and lon of RIGGS points (in deg,min,sec in .dat file)
RIGGSlat = -(RIGGS(:,4) + RIGGS(:,5)/60 + RIGGS(:,6)/(60*60));
RIGGSlon = RIGGS(:,7) + RIGGS(:,8)/60 + RIGGS(:,9)/(60*60);
RIGGSlon = - RIGGSlon .* RIGGS(:,10);  % RIGGS(:,10) is +1 if W, -1 if E

% add markers for RIGGS points
hold on, plot(RIGGSlon,RIGGSlat,'*k'), hold off

ccmap=get(gcf,'ColorMap'); ccmap(1,:)=[1 1 1]; ccmap(64,:)=[1 1 1];
set(gcf,'ColorMap',ccmap)
title('Color is speed in m/a.  Vectors shown are observed velocities at RIGGS points.')

