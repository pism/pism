function [U,dt] = diffusion(Lx,Ly,J,K,Dup,Ddown,Dright,Dleft,U0,tf)
% adaptive explicit method for non-constant diffusion
% equation  u_t = div (D grad u)
% form
%     U = diffusion(Lx,Ly,J,K,Dup,Ddown,Dright,Dleft,U0,tf)
% where Lx,Ly are half-widths of rectangular domain
%       D*    are (J-1) x (K-1) matrices for "staggered" grid
%       U0    is (J+1) x (K+1) matrix with initial values
% note: no error checking on sizes of D* or U0 matrices
% compare to heatadapt() for a first test:
%   >> J=50; K=50; D=ones(J-1,K-1);
%   >> [x,y]=ndgrid(-1:2/J:1,-1:2/K:1);
%   >> U0 = exp(-30*(x.*x + y.*y));
%   >> U = diffusion(1.0,1.0,J,K,D,D,D,D,U0,0.05);
%   >> surf(x,y,U);

% spatial grid and initial condition
dx = 2 * Lx / J;
dy = 2 * Ly / K;
[x,y] =  ndgrid(-Lx:dx:Lx, -Ly:dy:Ly); % (J+1) x (K+1) grid
U = U0;

% need  1 - 2 nu_x - 2 nu_y >= 0  :
maxD = max(max(Dup,Ddown),max(Dleft,Dright));
dt0 = 0.25 * min(dx,dy)^2 / max(max(maxD));
N = ceil(tf / dt0);
dt = tf / N;

mu_x = dt / (dx*dx);
mu_y = dt / (dy*dy);

% explicit steps
for n=1:N
   U(2:J,2:K) = U(2:J,2:K) + ...
       mu_y * Dup    .* ( U(2:J,3:K+1) - U(2:J,2:K) ) - ...
       mu_y * Ddown  .* ( U(2:J,2:K) - U(2:J,1:K-1) ) + ...
       mu_x * Dright .* ( U(3:J+1,2:K) - U(2:J,2:K) ) - ...
       mu_x * Dleft  .* ( U(2:J,2:K) - U(1:J-1,2:K) );
end

