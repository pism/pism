function exampletheta(J,h_level,lambdax,lambday)
% EXAMPLETHETA  Show the difference between "raw" and faster (Taylor P_4)
% methods for computing \theta in equation (49) in Schoof (2003).

if nargin < 1, J = 80; end
if nargin < 2, h_level = 1000.0; end
if nargin < 3, lambdax = 50e3; end
if nargin < 4, lambday = 50e3; end

L = 1200e3;
dx = 2 * L / J;
x = -L:dx:L;  y = x;  [xx,yy] = meshgrid(x, y);

% construct a bumpy bed
topg0 = 400 * sin(2 * pi * xx / 600e3) + 100 * sin(2 * pi * (xx + 1.5 * yy) / 40e3);

% OTHER BUMPYS:
%topg0 = 400 * sin(2 * pi * xx / 400e3) + 100 * randn(size(xx));
%topg0 = 100 * randn(size(xx));
%topg0 = zeros(size(xx));

[topg_smooth,theta,fasttheta] = gettheta(x,y,topg0,lambdax,lambday,h_level);

fprintf('min(topg_smooth) = %f, max(topg_smooth)     = %f\n',...
        min(min(topg_smooth)),max(max(topg_smooth)) )
fprintf('min(theta)       = %f, max(theta)     = %f\n',...
        min(min(theta)),max(max(theta)) )
fprintf('min(fasttheta)   = %f, max(fasttheta) = %f\n',...
        min(min(fasttheta)),max(max(fasttheta)) )

figure(1), mesh(x/1000,y/1000,topg0)
zmax = max(max(topg0));  zmin = min(min(topg0));
axis([-L/1000 L/1000 -L/1000 L/1000 zmin zmax])
xlabel('x  (km)'), ylabel('y  (km)'), zlabel('topg0  (m)')
title('initial topography')

figure(2), mesh(x/1000,y/1000,topg_smooth)
axis([-L/1000 L/1000 -L/1000 L/1000 zmin zmax])
xlabel('x  (km)'), ylabel('y  (km)'), zlabel('topg_smooth  (m)')
title('smoothed topography')

figure(3), clf
imagesc(x/1000,y/1000,theta,[min(min(theta)) max(max(theta))]), colorbar
xlabel('x  (km)'), ylabel('y  (km)')
title('theta (directly)')

figure(4), clf
imagesc(x/1000,y/1000,fasttheta,[min(min(theta)) max(max(theta))]), colorbar
xlabel('x  (km)'), ylabel('y  (km)')
title('approx theta (by fast method)')

difftheta = fasttheta-theta;
fprintf('max(difftheta) = %.6f, min(difftheta) = %.6f, average(difftheta) = %.6f\n',...
        max(max(difftheta)),min(min(difftheta)),mean(mean(difftheta)) )

figure(5), clf, imagesc(x/1000,y/1000,difftheta), colorbar
xlabel('x  (km)'), ylabel('y  (km)')
title('diff = fasttheta - theta')

