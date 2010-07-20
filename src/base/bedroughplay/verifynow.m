% VERIFYNOW  run a sequence of J=20,40,80,160 versions of tests B and L
% calls verifygeneralsia.m

disp('running Test B')
verifygeneralsia(20,'b')
verifygeneralsia(40,'b')
verifygeneralsia(80,'b')
verifygeneralsia(160,'b')

disp('running Test L (requires timestep control)')
verifygeneralsia(20,'l',4.0)
verifygeneralsia(40,'l',2.0)
verifygeneralsia(80,'l',1.0)
verifygeneralsia(160,'l',0.5)


%% SAMPLE OUTPUT FROM marmaduke.gi.alaska.edu, 11 July 2010:
%>> verifynow
%running Test B
%errors for 20 x 20 grid:
%average abs error in surface elev            = 22.310
%maximum abs error in surface elev            = 227.844
%errors for 40 x 40 grid:
%average abs error in surface elev            = 9.490
%maximum abs error in surface elev            = 241.464
%errors for 80 x 80 grid:
%average abs error in surface elev            = 2.801
%maximum abs error in surface elev            = 155.786
%errors for 160 x 160 grid:
%average abs error in surface elev            = 1.058
%maximum abs error in surface elev            = 109.437
%running Test L (requires timestep control)
%reading NetCDF file "testL_21.nc" to get exact solution for Test L ...
%errors for 20 x 20 grid:
%average abs error in surface elev            = 24.470
%maximum abs error in surface elev            = 359.084
%reading NetCDF file "testL_41.nc" to get exact solution for Test L ...
%errors for 40 x 40 grid:
%average abs error in surface elev            = 9.109
%maximum abs error in surface elev            = 359.084
%reading NetCDF file "testL_81.nc" to get exact solution for Test L ...
%errors for 80 x 80 grid:
%average abs error in surface elev            = 5.694
%maximum abs error in surface elev            = 318.349
%reading NetCDF file "testL_161.nc" to get exact solution for Test L ...
%errors for 160 x 160 grid:
%average abs error in surface elev            = 1.880
%maximum abs error in surface elev            = 235.668

