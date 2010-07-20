function H = halfar(t,x,y)
% compute the similarity solution, to the isothermal flat-bed SIA from Halfar (1983)
% constants as in Test B in Bueler et al (2005)
% inputs x,y can be matrices, of same size, in meters; scalar input t in seconds

H0 = 3600;    % constants in SI units
R0 = 750e3;
n = 3;
alpha = 1/9;
beta = 1/18;
secpera = 31556926;

Gamma = 9.0177e-13;  % m-3 s-1
t0 = (beta/Gamma) * (7/4)^3 * (R0^4/H0^7); % equation (9) in Bueler et al (2005); = 422.45 a

r = sqrt(x.*x + y.*y);
r = r / R0;  t = t / t0;
inside = max(0, 1 - (r / t^beta).^((n+1) / n));
H = H0 * inside.^(n / (2*n+1)) / t^alpha;

