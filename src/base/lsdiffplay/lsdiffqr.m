function lsdiffqr
% LSDIFFQR  Demo least squares by QR to compute a finite difference gradient 
% at the staggered grid point (x_i+1/2,y_j).  We have linear function of form
%   l(x,y) = p + q (x - x_i+1/2) + r (y - y_j)
% and we want it to go through 8 points  {(x,y,f(x,y))},  as shown in the 
% stencil figure lsdifffig.png.
%
% Having l(x,y) go through all 8 point is not possible, but we satisfy these
% 8 constraints approximately by least squares.  Once we find l(x,y) by least
% squares:     grad f = (q,r).

dx = 0.0001;      dy = 0.00014;       % sample values
hdx = 0.5*dx;   sdx = 1.5*dx;    % half and sesqui

% For next lines, refer to lsdifffig.png.
% Note that columns of A are already orthogonal (!) so QR decomposition
% is essentially trivial.  So we really can solve normal equations directly:
%   "A v = b" by   A' * W * A v = A' * W * b   where A'*W*A is diagonal (=R^2 below)
A = [1   -hdx   -dy;
     1    hdx   -dy;
     1   -sdx     0;
     1   -hdx     0;
     1    hdx     0;
     1    sdx     0;
     1   -hdx    dy;
     1    hdx    dy];
W = diag([1 1 1 1 1 1 1 1]);  % for now, equally weighting on all eight points
                              %   so W = eye(8)
b = [           f(-hdx,-dy) f(hdx,-dy) ...
     f(-sdx,0)  f(-hdx,0)   f(hdx,0)   f(sdx,0) ...
                f(-hdx,dy)  f(hdx,dy)]';

[Q,R] = qr(W*A,0);     % the so-called "reduced" QR where Q is 8 x 3 not 8 x 8
c = R \ (Q' * W * b);  % equivalently:  c = (A'*W*A) \ (A'*W*b)
%Q,R
%format rat, Q.^2, format long g

dfdx = c(2) % a bit different from smooth1D.m result
dfdy = c(3) % turns out to be Mahaffy

dfdxold = (b(5)-b(4))/dx % = centered
dfdyold = (b(7)+b(8)-b(1)-b(2))/(4*dy)  % = Mahaffy

  function z = f(x,y)  % sample function only
    z = sin(17 * (x-0.4)) * cos(20.3 * (y+1));
  end

