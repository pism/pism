function pism_matlab()
% creates a very boring dataset; an illustration of using MATLAB to create a PISM bootstrapping file.
% tested in MATLAB Version 7.7.0.471 (R2008b)

% set up the grid:
Lx = 1e6;
Ly = 1e6;
Mx = 101;
My = 201;
x = linspace(-Lx,Lx,Mx);
y = linspace(-Ly,Ly,My);

% create dummy fields
[xx,yy] = meshgrid(y,x);
topg = (xx + yy) / max(Lx, Ly);
acab = zeros(My,Mx);
artm = zeros(My,Mx) + 10.0; % 10 degrees Celsius
thk  = zeros(My,Mx) + 1.0; % 1 km thick


% create a file; NC_CLOBBER means "overwrite if exists"
ncid = netcdf.create('foo.nc', 'NC_CLOBBER');

% create dimensions:
x_id = netcdf.defDim(ncid, 'x', Mx);
y_id = netcdf.defDim(ncid, 'y', My);

% create coordinate variables:
x_var_id = netcdf.defVar(ncid, 'x', 'double', [x_id]);
y_var_id = netcdf.defVar(ncid, 'y', 'double', [y_id]);

% create variables corresponding to spatial fields:
acab_id = netcdf.defVar(ncid, 'acab', 'double', [x_id,y_id]);
artm_id= netcdf.defVar(ncid, 'artm', 'double', [x_id,y_id]);
topg_id = netcdf.defVar(ncid, 'topg', 'double', [x_id,y_id]);
thk_id = netcdf.defVar(ncid, 'thk', 'double', [x_id,y_id]);

% write attributes:
netcdf.putAtt(ncid, artm_id, 'units', 'Celsius');
netcdf.putAtt(ncid, thk_id, 'units', 'km');

% done defining dimensions and variables:
netcdf.endDef(ncid);

% write coordinate variables:
netcdf.putVar(ncid, x_var_id, x);
netcdf.putVar(ncid, y_var_id, y);

% write data itself:

% surface temperature:
netcdf.putVar(ncid, artm_id, artm);

% accumulation/ablation rate:
netcdf.putVar(ncid, acab_id, acab);

% bedrock elevation:
netcdf.putVar(ncid, topg_id, topg);

% ice thickness:
netcdf.putVar(ncid, thk_id, thk);

netcdf.close(ncid);

