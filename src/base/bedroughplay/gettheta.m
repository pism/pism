function [topgs,theta,fasttheta] = gettheta(x,y,topg0,lambdax,lambday,h_level)
% GETTHETA  Smooth the bed elevation map topg0, and then compute the correction
% factor "theta(D_0,X)" in equation (49) in Schoof (2003).  In PISM notation,
% theta = theta(h,x,y).  Recall that  Dnew = theta * D  is the corrected
% diffusivity.  The negative nth root of theta is given by this spatial average: 
%                       / lx / ly
% theta^{-1/n} = NORM * |    |  ( 1 - btilde(x,y,xi,nu) / (h - b0(x,y)) )^{-(n+2)/n} dxi dnu
%                       /-lx /-ly
% where
%   n = exponent in Glen law,
%   lx = lambdax, ly = lambday,
%   b0 = topg0,
%   h = h_level,   [a scalar here]
%   NORM = (4 * lambdax * lambday)^-1  [a normalizing factor; gives an average],
%   btilde(x,y,xi,nu) = bs(x+xi,y+nu) - b0(x,y),
% and the smoothed bed is
%                   / lx / ly
%  bs(x,y) = NORM * |    |     b0(x+xi,y+nu) dxi dnu
%                   /-lx /-ly
% with the same meaning for NORM.
% The output variables are
%   topgs = bs(x,y)       [smoothed bed; averages over (-lx,lx) x (-ly,ly) box],
%   theta = theta(h,x,y),
%   fasttheta = (approximation to theta using ??).


if ~all(size(topg0) == [length(x), length(y)]), error('array "topg" of wrong size'), end
if any(diff(diff(x))), error('  coordinate variable "x" not equally-spaced'), end
if any(diff(diff(y))), error('  coordinate variable "y" not equally-spaced'), end

n = 3.0; % Glen exponent

% constants needed below
pow = (n + 2) / n;
ccc2 = pow * (2 * n + 2) / (2 * n);
ccc3 = ccc2 * (3 * n + 2) / (3 * n);
ccc4 = ccc3 * (4 * n + 2) / (4 * n);
ccc5 = ccc4 * (5 * n + 2) / (5 * n);

J = length(x)-1;  K = length(y)-1;
dx = x(2) - x(1);  Nx = ceil(lambdax / dx)
dy = y(2) - y(1);  Ny = ceil(lambday / dy)

% smooth the bed; computation scales like  J * K * 2Nx * 2Ny
tic
topgs = zeros(J+1,K+1);
for j=1:J+1
  for k=1:K+1
    sm    = 0.0;
    count = 0;
    for r=-Nx:Nx
      for s=-Ny:Ny
        if (j+r >= 1) & (j+r <= J+1) & (k+s >= 1) & (k+s <= K+1)
          sm    = sm + topg0(j+r,k+s);
          count = count + 1;
        end
      end
    end
    topgs(j,k) = sm / count;
  end
end
fprintf('  [time to smooth the bed:                        %.5f s]\n',toc)

% get theta directly; again computation scales like  J * K * 2Nx * 2Ny
tic
theta = zeros(J+1,K+1); % note zero *is* the correct value for no valid pts in
                        % average below
for j=1:J+1
  for k=1:K+1
    topg_s = topgs(j,k);      % Schoof: "H"
    thk_s = h_level - topg_s; % Schoof: "D_0 - H"
    if thk_s > 0
      sm    = 0.0;
      maxtl = 0.0;
      count = 0;
      for r=-Nx:Nx
        for s=-Ny:Ny
          if (j+r >= 1) & (j+r <= J+1) & (k+s >= 1) & (k+s <= K+1)
            tl = topg0(j+r,k+s) - topg_s;  % Schoof: "h(x)"; has mean zero
                                           % over (j-Nx,j+Nx) x (k-Ny,k+Ny) patch
            maxtl = max(tl,maxtl);
            if thk_s > tl
              sm    = sm + ( 1 - tl / thk_s )^(-pow);
              count = count + 1;
            end
          end
        end
      end
      if thk_s > maxtl % Schoof: "theta is only defined if D_0-H is greater than the max of h"
        theta(j,k) = (sm / count)^(-n);
      end
    end
  end
end
fprintf('  [time to compute theta directly:                %.5f s]\n',toc)


% precompute coefficients in fasttheta method; this is just like smoothing, but
% we are averaging square and cube of local topography; again computation
% scales like  J * K * 2Nx * 2Ny
tic
maxtl = zeros(J+1,K+1);
C2 = zeros(J+1,K+1);
C3 = C2;
C4 = C2;
C5 = C2;
for j=1:J+1
  for k=1:K+1
    topg_s = topgs(j,k);      % Schoof: "H"
    maxtljk = 0.0;
    sm2   = 0.0;
    sm3   = 0.0;
    sm4   = 0.0;
    sm5   = 0.0;
    count = 0;
    for r=-Nx:Nx
      for s=-Ny:Ny
        if (j+r >= 1) & (j+r <= J+1) & (k+s >= 1) & (k+s <= K+1)
          tl    = topg0(j+r,k+s) - topg_s;  % Schoof: "h(x)"; has mean zero
                                            % over (j-Nx,j+Nx) x (k-Ny,k+Ny) patch
          tl2   = tl * tl;
          tl4   = tl2 * tl2;
          maxtljk = max(tl,maxtljk);
          sm2   = sm2 + tl2;
          sm3   = sm3 + tl2 * tl;
          sm4   = sm4 + tl4;
          sm5   = sm5 + tl4 * tl;
          count = count + 1;
        end
      end
    end
    C2(j,k) = sm2 / count;
    C3(j,k) = sm3 / count;
    C4(j,k) = sm4 / count;
    C5(j,k) = sm5 / count;
    maxtl(j,k) = maxtljk;
  end
end
C2 = ccc2 * C2;
C3 = ccc3 * C3;
C4 = ccc4 * C4;
C5 = ccc5 * C5;
fprintf('  [time to pre-compute coeffs maxtl,C2,C3,C4,C5:  %.5f s]\n',toc)


% get theta by fast method from stored maxtl,C2,C3 above; computation scales like  J * K
tic
fasttheta = zeros(J+1,K+1);
thks = h_level - topgs;
msk = (thks > maxtl);
thksm2 = zeros(J+1,K+1);  % just to allocate
thksm4 = zeros(J+1,K+1);  % just to allocate
thksm2(msk) = 1.0 ./ (thks(msk) .* thks(msk));
thksm4(msk) = thksm2(msk) .* thksm2(msk);
% high order version
fastomega = fasttheta;
fastomega(msk) = 1.0 + C2(msk) .* thksm2(msk) + C3(msk) .* (thksm2(msk) ./ thks(msk)) + ...
                         C4(msk) .* thksm4(msk) + C5(msk) .* (thksm4(msk) ./ thks(msk));
fastomega(fastomega < 0.001) = 0.001;
fasttheta(msk) = fastomega(msk).^(-3);
fasttheta(fasttheta > 1.0) = 1.0;
fasttheta(fasttheta < 0.0) = 0.0;
%lower order version:
%  fasttheta(msk) = ( 1.0 + C2(msk) .* thksm2(msk) + C3(msk) .* (thksm2(msk) ./ thks(msk)) ).^(-3);
fprintf('  [time to compute fasttheta from stored coeffs:  %.5f s]\n',toc)

