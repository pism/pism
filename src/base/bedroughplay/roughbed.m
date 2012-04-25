function roughbed(J,b_mag,b_wavelen,tfyears)
% ROUGHBED   Measure how the SIA code siageneral.m responds to rough bed topography.

if nargin < 1, J = 60; end
if nargin < 2, b_mag = 400; end
if nargin < 3, b_wavelen = 400e3; end
if nargin < 4, tfyears = 5000; end

% space-time grid parameters
L = 1200e3;
dx = 2 * L / J;
[x,y] = meshgrid(-L:dx:L, -L:dx:L);

% no mass balance
climatic_mass_balance = zeros(size(x));

% construct a bumpy bed
%topg = b_mag * randn(size(x));  % uncorrelated (unsmooth) random
topg = b_mag * sin(2 * pi * x / b_wavelen);  % pure sine

% initial surface over bumpy bed
secpera = 31556926;
t0 = 50;
thk0 = halfar(t0 * secpera,x,y);
h0 = thk0 + topg;
figure(1), mesh(x/1000,y/1000,h0)
zmax = max(max(h0));  zmin = min(min(h0));
axis([-L/1000 L/1000 -L/1000 L/1000 zmin zmax])
xlabel('x  (km)'), ylabel('y  (km)'), zlabel('h  (m)'), title('initial surface')

format long g
fprintf("  initial volume = %14.10e\n",sum(sum(thk0)) * dx * dx)

% run SIA code
dtyears = 10.0;
[h,dtlist] = siageneral(L,L,J,J,climatic_mass_balance,h0,topg,dtyears*secpera,tfyears*secpera);

thk = h - topg;
fprintf("  final volume   = %14.10e\n",sum(sum(thk)) * dx * dx)

% show final state
figure(2), mesh(x/1000,y/1000,h)
axis([-L/1000 L/1000 -L/1000 L/1000 zmin zmax])
xlabel('x  (km)'), ylabel('y  (km)'), zlabel('h  (m)'), title('final surface')

return

% adaptive time-stepping
figure(3),  plot(dtlist/secpera,'o')
xlabel('time steps  (count)'), ylabel('time-step dt  (a)')

