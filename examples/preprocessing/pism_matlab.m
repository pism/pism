function pism_matlab()
% PISM_MATLAB  Creates "from scratch" a boring dataset with the right format
% to use as a PISM bootstrapping file.  Example use of Matlab for this purpose.
%
% Usage, including a minimal PISM call to bootstrap from this file:
%    $ matlab
%    >> pism_matlab  % creates bar.nc
%    >> exit
%    $ pismr -i bar.nc -bootstrap -Mx 41 -My 41 -Mz 21 -Lz 4000 -Mbz 5 -Lbz 500 -y 1

% tested in MATLAB Version R2012b (8.0.0.783)

% set up the grid:
Lx = 1e6;
Ly = 1e6;
Mx = 51;
My = 71;
x = linspace(-Lx,Lx,Mx);
y = linspace(-Ly,Ly,My);

% create dummy fields
[xx,yy] = ndgrid(x,y);  % meshgrid() generates wrong ordering
acab = zeros(Mx,My);
artm = zeros(Mx,My) + 273.15 + 10.0; % 10 degrees Celsius
topg = 1000.0 + 200.0 * (xx + yy) / max(Lx, Ly);  % change "1000.0" to "0.0" to test
                                                  % flotation criterion, etc.
thk  = 3000.0 * (1.0 - 3.0 * (xx.^2 + yy.^2) / Lx^2);
thk(thk < 0.0) = 0.0;

% create a file; NC_CLOBBER means "overwrite if exists"
ncid = netcdf.create('bar.nc', 'NC_CLOBBER');

% create dimensions:
x_id = netcdf.defDim(ncid, 'x', Mx);
y_id = netcdf.defDim(ncid, 'y', My);

% create coordinate variables:
x_var_id = netcdf.defVar(ncid, 'x', 'float', [x_id]);
y_var_id = netcdf.defVar(ncid, 'y', 'float', [y_id]);

% write attributes:
netcdf.putAtt(ncid, x_var_id, 'long_name', 'easting');
netcdf.putAtt(ncid, x_var_id, 'standard_name', 'projection_x_coordinate');
netcdf.putAtt(ncid, y_var_id, 'long_name', 'northing');
netcdf.putAtt(ncid, y_var_id, 'standard_name', 'projection_y_coordinate');

% create variables corresponding to spatial fields:
% dimension transpose is standard: "float thk(y, x)" in NetCDF file
topg_id = netcdf.defVar(ncid, 'topg', 'float', [x_id,y_id]);
thk_id  = netcdf.defVar(ncid, 'thk',  'float', [x_id,y_id]);
acab_id = netcdf.defVar(ncid, 'climatic_mass_balance', 'float', [x_id,y_id]);
artm_id = netcdf.defVar(ncid, 'ice_surface_temp', 'float', [x_id,y_id]);

% write attributes:
netcdf.putAtt(ncid, topg_id, 'units', 'm');
netcdf.putAtt(ncid, topg_id, 'standard_name', 'bedrock_altitude');
netcdf.putAtt(ncid, thk_id, 'units', 'm');
netcdf.putAtt(ncid, thk_id, 'standard_name', 'land_ice_thickness');
netcdf.putAtt(ncid, acab_id, 'units', 'm year^-1');
netcdf.putAtt(ncid, acab_id, 'standard_name', 'land_ice_surface_specific_mass_balance');
netcdf.putAtt(ncid, artm_id, 'units', 'K');

% done defining dimensions and variables:
netcdf.endDef(ncid);

% write coordinate variables:
netcdf.putVar(ncid, x_var_id, x);
netcdf.putVar(ncid, y_var_id, y);

% surface temperature:
netcdf.putVar(ncid, artm_id, artm);

% accumulation/ablation rate:
netcdf.putVar(ncid, acab_id, acab);

% bedrock elevation:
netcdf.putVar(ncid, topg_id, topg);

% ice thickness:
netcdf.putVar(ncid, thk_id, thk);

netcdf.close(ncid);

disp('  PISM-bootable NetCDF file "bar.nc" written')
disp('  for example, run:')
disp('    $ pismr -i bar.nc -bootstrap -Mx 41 -My 41 -Mz 21 -Lz 4000 -Mbz 5 -Lbz 500 -y 1')

