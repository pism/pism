% ANTFIGS Plot various figures based on Antarctica runs of the ice model.
% Run "-of m" output from the ice model first.  (Example: ant_181_10yr.m.)
% Produces:
%    * map-plane view of temperature showing pressure melting mask;
%      assumes x,y,h,Tkd,pmMaskkd are defined.
%    * contour plot of semi-slice of temperature; assumes x,y,z,Tjd 
%      are defined

clear Mx dx hmask hmaskl
Mx=length(x);
My=length(y);
dx=(max(x)-min(x))/Mx;
hmask=double(h'>0);
zz=repmat(z',1,Mx);

Thomol = Tkd - (273.15 - H*8.66e-4);
Hid = repmat(H(:,(My-1)/2+1)',length(z),1);
Tidhomol = Tid - (273.15 - (Hid-zz)*8.66e-4);
Hjd = repmat(H((Mx-1)/2+1,:),length(z),1);
Tjdhomol = Tjd - (273.15 - (Hjd-zz)*8.66e-4);

% simple plot of basal homologous temp; white=pressure melting
hand1=figure;
imagesc(x,y,flipud(Thomol')), axis square, colorbar
M=get(hand1,'ColorMap');
M(64,:) = [1 1 1];
set(hand1,'ColorMap',M);

% % plot of basal homologous temp with superimpose pressure-melting mask
% figure
% surface(x-1.5*dx,y-1.5*dx,ones(Mx,Mx),0.0001*pmMaskkd',...
%     'AlphaData',pmMaskkd','FaceAlpha','flat','EdgeColor','none',...
%     'FaceColor',[0.4 0.4 0.4])
% hold on
% contourf('v6',x,y,Thomol',-35:5:0); % call if release 14 and later
% colorbar
% hold off
% axis tight, axis square, view(2)

% slice contour plot of temperature
figure
contourf(x,z,Tidhomol,-65:5:0);
colormap(jet), colorbar
hold on
fill([y(1) y y(end) y(1)],[z(end) h(:,(My-1)/2+1)' z(end) z(end)],'w')
title('temperature and profile semi-slice at x from -id');
xlabel('y (km)'); ylabel('z (m)'); 
hold off

% slice contour plot of temperature
figure
contourf(x,z,Tjdhomol,-55:5:0);
colormap(jet), colorbar
hold on
fill([x(1) x x(end) x(1)],[z(end) h((Mx-1)/2+1,:) z(end) z(end)],'w')
title('temperature and profile semi-slice at y from -jd');
xlabel('x (km)'); ylabel('z (m)'); 
hold off
