% SHOWTEMPATDEPTH illustrate formula used for put temp at depth

% constants
rhoi = 910.0; % kg m-3
c = 2009.0;   % J kg-1 K-1
ki = 2.10;    % W m-1 K-1
K = ki / (rhoi * c);  % m2 s-1

% data
Ts = 273.15 - 30;           % K      (-30 celsius)  surface temp
H =  2000;                  % m                     thickness
adot = 0.20 / 3.1556926e7;  % m s-1  (20 cm/year)   mass balance rate
g = 0.040;                  % W m-2  (40 mW/m^2)    geothermal flux

gamma = adot * H / K  % m s-1 m m-2 s = pure 
G = g * H / ki;        % W m-2 m W-1 m K = K

ccc = sqrt(gamma/2);
z = 0:H/200:H;
T = Ts + G * sqrt(pi/(2*gamma)) * (erf(ccc) - erf(ccc * z / H));
plot(T,z), grid on, xlabel T, ylabel z
