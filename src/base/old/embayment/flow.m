%function flow()
secpera = 3.1556926e7;
as = 1 / secpera;
rho = 910.0;
grav = 9.8;
%grav = 1/rho;
A = 4.0e-25;
B = A ^ (- 1 / 3);
Mx = 101;
My = 101;
xl = 0; xh = 500e3;
yl = -500e3; yh = 500e3;
dx = (xh - xl) / Mx;
dy = (yh - yl) / My;

[xx yy] = meshgrid(linspace(xl, xh, Mx), linspace(yl, yh, My));
x = xx / (xh - xl);
y = 2 * yy / (yh - yl);

uu = (0.5 * (1 + x.^2)) .* ((1.4 - y) .* (1.4 + y));
%uu = (x + 0.1) .* ((1.4 - y) .* (1.4 + y));
vv = - 0.4 * y - 0.6 * tan(y) .* (1 - x) - 0.5 * atan(y) .* (x + 0.5);
%vv = - 0.4 * y - 0.5 * tan(y) .* (1 - x) - 0.5 * atan(y) .* (x + 0.5);
uu = uu * 200 / secpera;
vv = vv * 200 / secpera;

theta  = -0.0 * pi / 2;
u = cos(theta) * uu - sin(theta) * vv;
v = sin(theta) * uu + cos(theta) * vv;

[u_x u_y] = gradient(u, dx, dy);
[v_x v_y] = gradient(v, dx, dy);

[u_xx u_xy] = gradient(u_x, dx, dy);

phi = (0.5 * u_x .^2 + 0.5 .* v_y .^ 2 + 0.5 .* (u_x + v_y) .^ 2 ...
         + 0.25 * (u_y + v_x) .^ 2);
[phi_x phi_y] = gradient(phi, dx, dy);
nu = (B / 2) .* phi .^ (-1 / 3);
[nu_x nu_y] = gradient(nu, dx, dy);

%H = 500 * (1 + 0.2 * (1 ./ (x + 2.5)) .* (1 + sin(1.0 .* y) .^ 2));
if (0)
  h = 500 * (0.8 + 0.3 * (1 ./ (x + 0.6)) ...
             .* (0.7 - 0.5*y.^2) + 0.4 * sin(1.0 .* y) .^ 2);
  b = 0 * x;
else
  %h = 500 * (3 ./ (x + 4) + 0.1 * abs(y) .^ 5 .* (1.2 - x));
  h = 500 * (3 ./ (x + 4) + 0.2 * abs(y) .^ 4 .* (1.2 - x));
  %h = 500 * (0.8 + 0.2 * (1 ./ (x + 0.6)) ...
             %.* (0.7 - 0.5*y.^2) + 0.3 * y .^ 2);
  b = -80 * x;
  %figure(2); contourf(h - b); colorbar; title('H');
  figure(3); contourf(nu); colorbar; title('nu');
end
H = h - b;
[h_x h_y] = gradient(h, dx, dy);

[f1, foo] = gradient(nu .* H .* (2 * u_x + v_y), dx, dy);
[g2, f2] = gradient(nu .* H .* (u_y + v_x), dx, dy);
[foo, g1] = gradient(nu .* H .* (u_x + 2 * v_y), dx, dy);
f = 2 * f1 + f2 - rho * grav * H .* h_x;
g = 2 * g1 + g2 - rho * grav * H .* h_y;

warning off MATLAB:divideByZero;
betaf = f ./ u;
betag = g ./ v;
betaf(abs(betaf) > 1e10) = 1e10;
betag(abs(betag) > 1e10) = 1e10;
warning on MATLAB:divideByZero;

if (0)
  ii = ceil(Mx / 2); jj = ceil(My / 2);
  B
  sprintf('u: (%f, %f) %e %e %e', x(ii,jj), y(ii,jj), secpera * u(ii,jj), secpera ...
          * u_x(ii,jj), secpera * u_xx(ii,jj))
  sprintf('phi: (%f, %f) %e %e %e', x(ii,jj), y(ii,jj), phi(ii,jj), ...
          phi_x(ii,jj), phi_y(ii,jj))
  sprintf('nu: (%f, %f) %e %e %e', x(ii,jj), y(ii,jj), nu(ii,jj), ...
          nu_x(ii,jj), nu_y(ii,jj))
  sprintf('rhsx: %e', f(ii,jj))
  sprintf('betax: %e', f(ii,jj) / u(ii,jj))
end

if (0)
  figure(1);
  subplot(2, 3, 1), contourf(xx, yy, u*secpera, 20), colorbar, title('u'), xlabel('x'), ylabel('y')
  subplot(2, 3, 2), contourf(xx, yy, v*secpera, 20), colorbar, title('v'), xlabel('x'), ylabel('y')
  %figure(3), surfl(xx, yy, h), title('Thickness'), xlabel('x'), ylabel('y')
  subplot(2, 3, 3), contourf(xx, yy, secpera * divergence(xx, yy, u .* H, v .* H)), colorbar, title('div(flux)'), xlabel('x'), ylabel('y')

  subplot(2, 3, 4), contourf(xx, yy, betaf, 15), colorbar, title('beta (x)'), xlabel('x'), ylabel('y')
  subplot(2, 3, 5), contourf(xx, yy, betag, 15), colorbar, title('beta (y)'), xlabel('x'), ylabel('y')
  subplot(2, 3, 6), contourf(xx, yy, betaf - betag, 15), ...
      colorbar, title('beta (y) - beta (x)'), xlabel('x'), ylabel('y')

  figure(2);
  subplot(1, 3, 1), contourf(xx, yy, f, 15), ...
      colorbar, title('f'), xlabel('x'), ylabel('y')
  subplot(1, 3, 2), contourf(xx, yy, g, 15), ...
      colorbar, title('g'), xlabel('x'), ylabel('y')
  subplot(1, 3, 3), contourf(xx, yy, f - g, 15), ...
      colorbar, title('f - g'), xlabel('x'), ylabel('y')
end

if (1)
  figure(2); clf; axis([0, Mx, -My, My]); axis equal; s = 8;
  subplot(1,2,1); contourf(xx/1e3, yy/1e3, secpera*sqrt(u.^2 + v.^2), 40); shading flat; view(2);
  hold on;
  quiver(xx(1:s:Mx, 1:s:My)/1e3, ...
         yy(1:s:Mx, 1:s:My)/1e3, ...
         secpera*u(1:s:Mx, 1:s:My), ...
         secpera*v(1:s:Mx, 1:s:My), 0.8);
  colorbar; title('Velocity'); xlabel('x'), ylabel('y');
  subplot(1,2,2); contourf(xx/1e3, yy/1e3, H, 40); shading flat; view(2);
  hold on; quiver(xx(1:s:Mx, 1:s:My)/1e3, ...
                  yy(1:s:Mx, 1:s:My)/1e3, ...
                  -h_x(1:s:Mx, 1:s:My), ...
                  -h_y(1:s:Mx, 1:s:My));
  colorbar; title('Thickness'), xlabel('x'), ylabel('y')
end


if (1)
  figure(1); clf; axis([0, Mx, -My, My]); axis equal; s = 8;
  beta1 = betaf;
  beta2 = betag; beta2(betag ~= betag) = 0;
  subplot(1,2,1); contourf(xx/1e3, yy/1e3, beta1, 40); shading flat; view(2);
  colorbar; title('\beta_1'); xlabel('x'), ylabel('y');
  subplot(1,2,2); contourf(xx/1e3, yy/1e3, beta2, 40); shading flat; view(2);
  colorbar; title('\beta_2'); xlabel('x'), ylabel('y');
end
