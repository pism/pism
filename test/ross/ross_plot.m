% ROSS_PLOT Plots computed speed of Ross and shows observed values from RIGGS
% data.  
%    Assumes foo.m, the result of PISM Ross computation, has already been
% run (i.e. "$ pisms -ross -o foo -of m" and then do ">> foo"), so variables 
% c,ubar,vbar,H,mask,x,y are defined.  Reads RIGGS data from riggs_ELBclean.dat.
% Compare figures 3 from (Bently 1984) and figures 11 and 14 from 
% (Thomas et al 1984) and figures in (MacAyeal et al 1996).
%    Reports final values 'ChiSqr' and 'max_computed_speed' which can (I think)
% be compared to Table 1 in (MacAyeal et al 1996).
% ELB 2/4/07; 2/11/07

% see 111by147.dat for these ranges
dlat = (-5.42445 - (-12.3325))/110;
gridlatext = linspace(-12.3325 - dlat * 46,-5.42445,147);
gridlon = linspace(-5.26168,3.72207,147);

% load RIGGS data FROM D. MACAYEAL TO ELB ON 19 DEC 2006.
load -ascii riggs_ELBclean.dat
RIGGS=riggs_ELBclean;
clear riggs_ELBclean;

% show computed speed as color
cforplot=c;  cforplot(H<20) = -20; cforplot(mask==1) = -20;
figure
imagesc(gridlon,gridlatext,cforplot'), colorbar
h=get(gcf,'CurrentAxes');  set(h, 'YDir', 'normal')
axis([-5.26168 3.72207 -13 -5.42445])
xlabel('RIGGS grid longitude (deg E)'), ylabel('RIGGS grid latitude (deg N)')
ccmap=get(gcf,'ColorMap');  ccmap(1,:)=[1 1 1];  set(gcf,'ColorMap',ccmap)

% compute grid lat and lon of RIGGS points (in deg,min,sec in .dat file); 
% throw out the ones which are not in model domain; 132 remain
RIGGSlat = -(RIGGS(:,4) + RIGGS(:,5)/60 + RIGGS(:,6)/(60*60));
RIGGSlon = RIGGS(:,7) + RIGGS(:,8)/60 + RIGGS(:,9)/(60*60);
RIGGSlon = - RIGGSlon .* RIGGS(:,10);  % RIGGS(:,10) is +1 if W, -1 if E
cRIGGS = griddata(gridlon,gridlatext,cforplot',RIGGSlon,RIGGSlat,'nearest');
rig = RIGGS(cRIGGS>0,:); 
riglon = RIGGSlon(cRIGGS>0); riglat = RIGGSlat(cRIGGS>0); 

% add markers for RIGGS points, then quiver observed velocities
hold on, plot(riglon,riglat,'*k')
rigu = sin((pi/180)*rig(:,13)) .* rig(:,11);
rigv = cos((pi/180)*rig(:,13)) .* rig(:,11);
quiver(riglon,riglat,rigu,rigv,1,'k');

% quiver the computed velocities at the same points; note reversal of u,v in model
spera = 31556926;
uATrig = spera * griddata(gridlon,gridlatext,vbar',riglon,riglat,'linear');
vATrig = spera * griddata(gridlon,gridlatext,ubar',riglon,riglat,'linear');
quiver(riglon,riglat,uATrig,vATrig,max(max(c))/max(max(sqrt(rigu.^2 + rigv.^2))),'r');
%quiver([riglon riglon],[riglat riglat],[rigu uATrig],[rigv vATrig],'k');
hold off
title('Color is speed in m/a.  Arrows are observed (black) and computed (red) velocities at RIGGS points.')

% report results comparable to Table 1 in (MacAyeal et al 1996)
ChiSqrActual = sum( ((uATrig - rigu).^2 + (vATrig - rigv).^2) / (30^2) );
ChiSqr = ChiSqrActual * (156/132)
max_computed_speed = max(max(c))

% show observed versus computed scatter plot as in Figure 2 in (MacAyeal et al 1996)
figure
plot(sqrt(uATrig.^2 + vATrig.^2),sqrt(rigu.^2 + rigv.^2),'.k','Markersize',12)
hold on, plot([0 1000],[0 1000],'LineWidth',2), hold off
