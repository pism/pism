% ICEMODELFIGS Plot various figures used in 
%    Bueler, Kallen-Brown,  and Lingle (2006) "Exact solutions to 
%    thermocoupled ice-sheet models: a new standard for verification," 
%    to be submitted to J. Glaciol.
% Run "-of m" output from the ice model first.  (Example: eisF_121_201.m.)
% Produces:
%    * map-plane view of temperature showing pressure melting mask as in 
%      (Saito et al 2006); assumes x,y,h,Tkd,pmMaskkd are defined.
%    * map-plane view of L = |dSigma/dT|, with marginal junk removed;
%      assumes x,y,h,Sigmakd are defined.
%    * contour plot of semi-slice of temperature; assumes x,y,z,Tjd 
%      are defined

figure
clear Mx dx hmask hmaskl
Mx=length(x);
dx=(max(x)-min(x))/Mx;
hmask=double(h>0);

% old version of contour DOES work, new not!
%if (version('-release') <= 13)
    [c,hand]=contour(x,y,Tkd,200:2:300,'k');    
%else
%    [c,hand]=contour('v6',x,y,Tkd,200:2:300,'k');
%end
%clabel(c,hand)
colormap gray
surface(x-1.5*dx,y-1.5*dx,ones(Mx,Mx),32*pmMaskkd,...
    'AlphaData',pmMaskkd,'FaceAlpha','flat','EdgeColor','none',...
    'FaceColor',[0.4 0.4 0.4])
%    'CDataMapping','direct')
xl=[x x(end)+dx x(end)+2*dx]; yl=[y y(end)+dx y(end)+2*dx];
hmaskl=zeros(Mx+2,Mx+2);
hmaskl(1:Mx,1:Mx)=hmask;
surface(xl-1.5*dx,yl-1.5*dx,ones(Mx+2,Mx+2),64*(1-hmaskl),...
    'AlphaData',(1-hmaskl),'FaceAlpha','flat','EdgeColor','none',...
    'FaceColor','white')
axis square, view(2)

% next part adds a figure showing Lip = |dSigma/dT| for the purpose of
% analyzing the spokes
R=8.314; % J/(mol K); (EISMINT II value in Paterson-Budd)
Qcold=6.0e4; % J/mol; (EISMINT II value in Paterson-Budd for T<263)
Qwarm=13.9e4; % J/mol; (EISMINT II value in Paterson-Budd for T>263)
Sigm=Sigmakd;
% note following rule: If SigmaCkd is defined, i.e. the compensatory Sigma,
% then use simple Arrhenius.  Otherwise use Paterson-Budd.  Removes
% marginal junk by slightly different rules in the two cases.
if exist('SigmaCkd')
    Sigm( h<1000 ) = 0;
    Lip = (Qcold/R)*abs(Sigm).*Tkd.^(-2);
else
    Sigm( h<1200 ) = 0;
    Lip = (Qcold/R)*abs(Sigm).*Tkd.^(-2);
    Lip(Tkd >= 263) = (Qwarm/Qcold)*Lip(Tkd >= 263);
end
figure
%if (version('-release') <= 13)
    contour(x,y,Lip)
%else
%    contour('v6',x,y,Lip)
%end
surface(x,y,Lip,'LineStyle','none',...
   'AlphaData',ones(size(Sigm)))
colormap gray, axis square, view(2)

return
% slice contour plot of temperature
figure
[cc chand]=contour(x((Mx-1)/2+1:Mx),z,Tjd(:,(Mx-1)/2+1:Mx),...
              [220 230 240 250 260 270 280],'k-');
hold on, clabel(cc,chand,'FontSize',14)
[cc2 chand2]=contour(x((Mx-1)/2+1:Mx),z,Tjd(:,(Mx-1)/2+1:Mx),...
              224:2:228,'k:','LineWidth',0.5);
clabel(cc2,chand2,'FontSize',12)
fill([x((Mx-1)/2+1) x((Mx-1)/2+1:Mx) x(end) x((Mx-1)/2+1)],...
    [z(end) h((Mx-1)/2+1:Mx,(Mx-1)/2+1)' z(end) z(end)],'w')
title('temperature and profile semi-slice at y from -jd');
xlabel('x (km)'); ylabel('z (m)'); 
axis([0 (8/9)*x(end) 0 z(end)])
hold off
