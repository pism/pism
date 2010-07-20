% TESTPADE  Compare Pade and MacLaurin expansions of a function  f(z)  similar
% to that which appears in Schoof (2003) for the definition of \theta.

k = 5/3;
alpha = 1.5;  % Pade looks good
%alpha = 1.05;   % Pade generates new pole!!

z = -1:0.01:1;
f = (1 - z/alpha).^(-k);

a1 = k / alpha;
a2 = a1 * (k + 1) / (2 * alpha);
a3 = a2 * (k + 2) / (3 * alpha);
a4 = a3 * (k + 3) / (4 * alpha);
a5 = a4 * (k + 4) / (5 * alpha);

mac3 = 1 + a1 * z + a2 * z.^2 + a3 * z.^3;

mac5 = mac3 + a4 * z.^4 + a5 * z.^5;

% these come from system of equations:
%    p1      - a0 q1 = a1
%         p2 - a1 q1 = a2
%            - a2 q1 = a3
q1 = - (k + 2) / (3 * alpha);  % = - a3 / a2
p1 = a1 + q1;
p2 = a2 + a1 * q1;

pad = (1 + p1 * z + p2 * z.^2) ./ (1 + q1 * z);

plot(z,f,z,mac3,z,mac5,z,pad)
axis([-1, 1, 0.5*min(f), 1.3*max(f)])
legend('f(z)','maclaurin O(3)','maclaurin O(5)','pade')

