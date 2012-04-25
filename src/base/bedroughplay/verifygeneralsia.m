function verifygeneralsia(J,test,dtyears)
% VERIFYGENERALSIA   Verify siageneral() using exact solutions.  Test naming
% is from PISM.
%
% call for Test B:  verifygeneralsia(J)
% Uses the Halfar (1983) similarity solution.  Runs from t=200a to t=20000a.
% Uses "outer" time steps of 10.0 years.
% where  J       = (number of grid spaces in x,y directions)
%
% call for Test L:  verifygeneralsia(J,'L')
% where  J       = (number of grid spaces in x,y directions)
% Here J must be from list {20,40,80,160}.  Run generate_testL.sh first; uses
% stored values of the exact solution.  Uses "outer" time steps of 10 years.
%
% general call:   verifygeneralsia(J,test,dtyears)
% where  J       = (number of grid spaces in x,y directions); list above if test=='L'
%        test    = one of {'B','L'}; if lower case then suppress dots
%        dtyears = "outer" timestep in years; diffusion.m can cut into smaller

if nargin<1, J=40; end  % default grid spaces

if nargin<2, doHalfar = true; dodots = true;
else
  if test == 'B', doHalfar = true; dodots = true;
  elseif test == 'b', doHalfar = true; dodots = false;
  elseif test == 'L', doHalfar = false; dodots = true;
  elseif test == 'l', doHalfar = false; dodots = false;
  else, error('argument "test" must be B or L'), end
end

if nargin<3, dtyears = 10.0; end;

if doHalfar
  L = 1200e3;
else
  L = 900e3;
  if ~any(J == [20 40 80 160]), error('for Test L, J must be one of {20,40,80,160}'), end
end

secpera = 31556926.0;
dx = 2 * L / J;
[x,y] = meshgrid(-L:dx:L, -L:dx:L);

if doHalfar
  t0 = 200.0;
  t1 = 20000.0;
  h0 = halfar(t0 * secpera,x,y);
  h1exact = halfar(t1 * secpera,x,y);
  topg = zeros(size(h0));  % flat bed
  a = topg;                % no surface balance
else
  filename = ['testL_' num2str(J+1) '.nc'];
  fprintf('reading NetCDF file "%s" to get exact solution for Test L ...\n', filename)
  nc_id = netcdf.open(filename, 'NOWRITE' );  % ONLY WORKS IN MATLAB
  climatic_mass_balance_id = netcdf.inqVarID(nc_id, 'climatic_mass_balance');
  topg_id = netcdf.inqVarID(nc_id, 'topg');
  usurf_id = netcdf.inqVarID(nc_id, 'usurf');
  a = netcdf.getVar(nc_id, climatic_mass_balance_id, 'double') / secpera;
  topg = netcdf.getVar(nc_id, topg_id, 'double');
  h0 = netcdf.getVar(nc_id, usurf_id, 'double');
  netcdf.close(nc_id);
  t0 = 0.0;
  t1 = 5000.0;
  h1exact = h0;
end

runtime = (t1 - t0) * secpera;
if ~dodots, runtime = -runtime; end % sign of this arg to siageneral() detemines dots

[h1approx,dtlist] = siageneral(L,L,J,J,a,h0,topg,dtyears * secpera,runtime);

err = h1approx - h1exact;
fprintf('errors for %d x %d grid:\n',J,J)
fprintf('average abs error in surface elev            = %.3f\n',...
        sum(sum(abs(err)))/(J+1)^2 )
fprintf('maximum abs error in surface elev            = %.3f\n',...
        max(max(abs(err))) )

if ~dodots, return, end

figure(1), imagesc((-L:dx:L)/1000.0, (-L:dx:L)/1000.0, err)
shading('flat'), axis square, colorbar, xlabel('x  (km)'), ylabel('y  (km)')
title('surface elevation error  (m)')

% figure showing adaptive time-stepping:
figure(3), plot(dtlist / secpera,'o')
xlabel('step'), ylabel('length of step in years')

% side-by-side comparison of numerical and exact result:
figure(2)
subplot(121), surf(x/1000,y/1000,h1exact), shading('flat')
xlabel('x (km)'), ylabel('y (km)'), zlabel('exact surface elevation (m)')
subplot(122), surf(x/1000,y/1000,h1approx), shading('flat')
xlabel('x (km)'), ylabel('y (km)'), zlabel('numerical surface elevation (m)')

