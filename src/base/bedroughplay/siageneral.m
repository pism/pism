function [h,dtlist] = siageneral(Lx,Ly,J,K,a,h0,b,deltat,tf)
% SIAGENERAL   Numerical solution of isothermal n=3 SIA with non-flat bed
% and general time-independent mass balance.  Uses Mahaffy (1976) method to
% evaluate map-plane diffusivity.  Calls diffusion.m, which does adaptive
% explicit time-stepping within the major time step.
%
% equation:   h_t = a + div (D grad h)
% where
%   a     = surface mass balance     [ m s^-1    ]
%   h     = ice surface elevation    [ m         ]
%   b     = bedrock elevation        [ m         ]
%   D     = Gamma (h-b)^5 |grad h|^2 [ m^2 s^-1  ]
%   Gamma = 2/5 A (rho g)^3          [ m^-3 s^-1 ]
%
% call:
%   [h,dtlist] = siageneral(Lx,Ly,J,K,a,h0,b,deltat,tf)
%
% inputs:
%   Lx,Lx  = half lengths of rectangle in x,y directions
%   J,K    = number of subintervals in x,y directions
%   a      = surface mass balance, a (J+1)x(K+1) array
%   h0     = initial surface elevation, a (J+1)x(K+1) array
%   b      = bedrock elevation, a (J+1)x(K+1) array
%   deltat = major time step
%   tf     = final time where runs start from t=0; if negative then suppress dots at stdout
%
% outputs:
%   h      = numerical approx of surface elevation at final time
%   dtlist = list of time steps used adaptively in diffusion.m
% 
% example: see verifygeneralsia.m

% constants
g = 9.81;
rho = 910.0;
secpera = 31556926;
A = 1.0e-16 / secpera;
Gamma  = 2 * A * (rho * g)^3 / 5; % see Bueler et al (2005)
h = h0;

if tf < 0.0, tf = -tf; dodots = false;
else, dodots = true; end

dx = 2 * Lx / J;  dy = 2 * Ly / K;
N = ceil(tf / deltat);  deltat = tf / N;
% indexing:
j  = 2:J;    k = 2:K;
nk = 3:K+1; sk = 1:K-1; ej = 3:J+1; wj = 1:J-1; % north,south,east,west

t = 0;
dtlist = [];
for n=1:N
  % regular grid thicknesses
  HH = h - b;
  HH(h < b) = 0.0;
  % compute staggered grid thicknesses
  HHup  = 0.5 * ( HH(j,nk) + HH(j,k) ); % up
  HHdn  = 0.5 * ( HH(j,k) + HH(j,sk) ); % down
  HHrt  = 0.5 * ( HH(ej,k) + HH(j,k) ); % right
  HHlt  = 0.5 * ( HH(j,k) + HH(wj,k) ); % left
  % staggered grid value of |grad h|^2 = (slope squared)
  a2up = (h(ej,nk) + h(ej,k) - h(wj,nk) - h(wj,k)).^2 / (4*dx)^2 + ...
         (h(j,nk) - h(j,k)).^2 / dy^2;
  a2dn = (h(ej,k) + h(ej,sk) - h(wj,k) - h(wj,sk)).^2 / (4*dx)^2 + ...
         (h(j,k) - h(j,sk)).^2 / dy^2;
  a2rt = (h(ej,k) - h(j,k)).^2 / dx^2 + ...
         (h(ej,nk) + h(j,nk) - h(ej,sk) - h(j,sk)).^2 / (4*dy)^2;
  a2lt = (h(j,k) - h(wj,k)).^2 / dx^2 + ...
         (h(wj,nk) + h(j,nk) - h(wj,sk) - h(j,sk)).^2 / (4*dy)^2;     
  % Mahaffy evaluation of staggered grid diffusivity D
  Dup  = Gamma * HHup.^5 .* a2up;
  Ddn  = Gamma * HHdn.^5 .* a2dn;
  Drt  = Gamma * HHrt.^5 .* a2rt;
  Dlt  = Gamma * HHlt.^5 .* a2lt;
  % call *adaptive* diffusion() to do time step; returns its time step
  [h,dtadapt] = diffusion(Lx,Ly,J,K,Dup,Ddn,Drt,Dlt,h,deltat);
  % PDE splitting: now add mass balance
  h = h + deltat * a;
  % now enforce b.c.
  h(h < b) = b(h < b);
  if dodots, fprintf('.'), end
  t = t + deltat;
  dtlist = [dtlist dtadapt];
end
if dodots, fprintf('\n'), end

% plot the diffusivity for the final state:
%x = linspace(-Lx+dx,Lx-dx,J-1);  y = linspace(-Ly+dy,Ly-dy,K-1);
%figure(99), surf(x,y,0.25*(Dup+Ddn+Drt+Dlt))

