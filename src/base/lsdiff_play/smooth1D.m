function [x, f, df_new, df_old, xs, df_news, df_olds] = smooth1D(N,fcn)
% SMOOTH1D   Experiment with 5 pt stencil for reg, and 4 pt for staggered,
% differentiation, but in only one dimension.  Using f(x) = sin(4 pi x) on
% 0 <= x <= 1 and periodic grid.  Comparing old CENTERED and new LEAST SQUARES
% methods.
%
% The new method works this way for regular grid:  At x=x_j we have line of
% form
%      l(x) = p + q (x - x_j).
% We "fit" this line to the five-point "data"
%      ( x_j-2, f(x_j-2) )
%      ( x_j-1, f(x_j-1) )
%      ( x_j,   f(x_j)   )
%      ( x_j+1, f(x_j+1) )
%      ( x_j+2, f(x_j+2) )
% Then we find  l'(x_j) = q.
% The "fitting" is by the standard QR method for solving least squares problems.
% A non-trivial weighting (diagonal matrix W) of the five points is allowed, but
% for now they are all equally-weighted.
%
% The new staggered method uses the four points at  x_j-1, x_j, x_j+1, x_j+2
% to find f'(x_j+1/2).  Also l(x) = p + q (x - x_j+1/2).  Otherwise the method
% is the same.
%
% usage:     [x, f, df_new, df_old, xs, df_news, df_olds] = smooth1D(N)
%
% example 1: case of N large enough to have some chance of extracting slopes by 
% both old and new methods
%   >> N=20; [x, f, new, old, xs, news, olds] = smooth1D(N);
%   >> plot(x,10*f,x,new,x,old,xs,news,xs,olds)
%   >> legend('10 f(x)','new','old','new stag','old stag')
%
% example 2: now try with N=9 to see how the new method succeeds in damping
% the derivatives in a case which is "hopeless" by centered-differencing
%
% example 3: now try with N=40 to see how the new and old are same for smooth

if N < 5, error('N >= 5 required'), end

dx = 1.0/N;
x  = 0:dx:1-dx;  % periodic grid
xs = x+(dx/2);

if nargin < 2,  f = sin(4*pi*x);
else,           f = feval(fcn,x);  end

% for NEW REGULAR:
% we minimize ||W (A v - b)||_2  by solving  W A v = W b  by least squares by QR
% in more serious version, matrices inv(R) and (Q' * W) would be pre-stored 
W = diag([1 1 1 1 1]);  % for now, all points count equally in 5 pt stencil
A = [1, -2*dx;  % "Av = b" would mean straight line l(x) = p + q (x - x_j)
     1,   -dx;  % would go through all five points
     1,     0;
     1,    dx;
     1,  2*dx];
[Q,R] = qr(W*A,0);  % so-called "reduced" QR where Q is 5 x 2 not 5 x 5

% for NEW STAGGERED:
Ws = diag([1 1 1 1]);
As =[1, -1.5*dx;
     1, -0.5*dx;
     1,  0.5*dx;
     1,  1.5*dx];
[Qs,Rs] = qr(Ws*As,0);

df_old=x; df_new=x; df_olds=xs; df_news=xs; % just to allocate

for j=1:N
  % REGULAR: OLD way is centered-differencing
  if j==1
    df_old(j) = (f(j+1) -  f(N) ) / (2*dx);
  elseif j==N
    df_old(j) = ( f(1)  - f(j-1)) / (2*dx);
  else % general case
    df_old(j) = (f(j+1) - f(j-1)) / (2*dx);
  end

  % REGULAR: NEW way is least-squares
  if j==1
    b = [f(N-1)  f(N)  f(j) f(j+1) f(j+2)]';
  elseif j==2
    b = [ f(N)  f(j-1) f(j) f(j+1) f(j+2)]';
  elseif j==N-1
    b = [f(j-2) f(j-1) f(j) f(j+1)  f(1) ]';
  elseif j==N
    b = [f(j-2) f(j-1) f(j)  f(1)   f(2) ]';
  else % general case
    b = [f(j-2) f(j-1) f(j) f(j+1) f(j+2)]';
  end
  c = R \ (Q' * W * b);
  df_new(j) = c(2);  % l(x) = c(1) + c(2) * (x - x_j)  so f'(x_j) ~ c(2)

  % STAGGERED: OLD way is centered-differencing
  if j==N
    df_olds(j) = ( f(1)  - f(j)) / dx;
  else % general case
    df_olds(j) = (f(j+1) - f(j)) / dx;
  end

  % STAGGERED: NEW way is least-squares
  if j==1
    bs = [ f(N)  f(j) f(j+1) f(j+2)]';
  elseif j==N-1
    bs = [f(j-1) f(j) f(j+1)  f(1) ]';
  elseif j==N
    bs = [f(j-1) f(j)  f(1)   f(2) ]';
  else % general case
    bs = [f(j-1) f(j) f(j+1) f(j+2)]';
  end
  cs = Rs \ (Qs' * Ws * bs);
  df_news(j) = cs(2);  % l(x) = cs(1) + cs(2) * (x - x_j+1/2)  so f'(x_j+1/2) ~ cs(2)

end

