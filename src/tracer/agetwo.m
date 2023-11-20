function [x,a1,a2,v,t,a1sum,a2sum] = agetwo(N,Tf)
% AGETWO Compare two methods for solving a pure advection problem, which
% can be thought of as the "conservative form" age equation
%   (1)    a_t + (v a)_x = 1
% on the interval 0 <= x <= 10.  Here  a = a(x,t)  is the age and  v = v(x)
% is the scalar velocity.
%
% Note that in 3D ice problems we are solving  a_t + div((u,v,w) a) = 1
% which is equivalent to  a_t + (u,v,w) . grad(a) = 1  because the velocity
% vector field (u,v,w) is incompressible.  In this example code we are
% paying attention to conservative properties of tracer schemes.  The age
% interpretation is incidental.
%
% Both schemes are upwinded explicit methods.  The first scheme applies to
% the equivalent equation
%   (2)    a_t + v a_x = 1 - v' a
% and uses classical upwinding for the left side.  The second scheme uses
% the conservative form (1) stated above.  Both schemes use CFL to determine
% the time step.
%
% The particular example problem is on the interval  [0,L] = [0,10]  and has
%    v(x) = cos(pi x / (2L)) (1 - 0.7 exp(-(x-3)^2))
% which is positive but has  v(L) = 0.  Roughly, the left end  x = 0  is
% analogous to the ice surface, and the right end  x = L  is analogous to
% the ice base, and  v(L) = 0  is analogous to a frozen base.  Thus
%   PDE   equation (1) above
%   BC    a(0,t) = 0
%   IC    a(x,0) = 0
%
% Try:  COMPAREAGE

L = 10;
dx = L/N;
x = 0:dx:L;
v = vv(x,L);
a1 = zeros(size(x));  % zero initial condition
a2 = a1;

dt = dx / max(vv(x,L));
M = ceil(Tf / dt);
dt = Tf / M;
fprintf('  M=%d and dt=%.4f\n',M,dt)
t = 0:dt:Tf;  % length M+1
a1sum = zeros(size(t));  a2sum = a1sum;

nu = dt / dx;

for m = 1:M
  % method one: classic first-order upwinding
  a1new(1) = 0.0;  % left b.c.  a(0,t) = 0
  for j = 2:N+1
    vj = vv(x(j),L);
    if vj > 0
      da = a1(j) - a1(j-1);
    else
      da = a1(j+1) - a1(j);
    end
    a1new(j) = a1(j) - nu * vj * da + dt - dt * dvv(x(j),L) * a1(j);
  end
  a1 = a1new;
  a1sum(m+1) = (dx/2) * [1 repmat([2],1,N-1) 1] * a1';  % trap rule (finite vol)

  % method two: conservative first-order upwinding
  % precompute cell-boundary fluxes  q = v a
  q = zeros(1,N);
  for j = 1:N
    vmid = 0.5 * (vv(x(j),L) + vv(x(j+1),L));
    if vmid > 0
      q(j) = vmid * a2(j);
    else
      q(j) = vmid * a2(j+1);
    end
  end
  % update a2
  % left b.c.  a(0,t) = 0 so q(0,t) = 0
  % q(x=0) = 0 but there *is* a dx/2 cell:
  a2new(1) = a2(1) - 2 * nu * (q(1) - 0) + dt;
  for j = 2:N
    a2new(j) = a2(j) - nu * (q(j) - q(j-1)) + dt;
  end
  % need q(x=L); there *is* a dx/2 cell:
  a2new(N+1) = a2(N+1) - 2 * nu * (vv(x(N+1),L) * a2(N+1) - q(N)) + dt;
  a2 = a2new;
  a2sum(m+1) = (dx/2) * [1 repmat([2],1,N-1) 1] * a2';  % trap rule (finite vol)
end

  function y = vv(x,L)
  y = cos(pi * x / (2*L)) .* (1 - 0.7 * exp(-(x-3).^2));
  end

  function dy = dvv(x,L)
  ee = exp(-(x-3).^2);
  cc = cos(pi * x / (2*L));
  ss = sin(pi * x / (2*L));
  dy = - (pi/(2*L)) * ss .* (1 - 0.7 * ee) - 0.7 * cc .* (-2 * (x-3)) .* ee;
  end

end
