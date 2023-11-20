function sumerrs = compareage(N)
% COMPAREAGE  runs agetwo.m with goal of comparing first-order upwinding

L = 10;
Tf = 10.0;
[x,a1,a2,v,t,a1sum,a2sum] = agetwo(N,Tf);
sumerrs = [a1sum(end)-L*Tf; a2sum(end)-L*Tf];

figure(1), plot(x,v), grid on
title('scalar velocity v(x)')
xlabel x, ylabel v

set(0,'defaultlinemarkersize',8)
figure(2), plot(x,a1,'-o',x,a2,'-o'), grid on
legend('first method','second method (conserving)')
title(sprintf('solution a(x,t) at t=%.4f',Tf))
xlabel x, ylabel('age  a(x,t)')

figure(3), plot(t,a1sum-L*t,t,a2sum-L*t), grid on
legend('first method','second method (conserving)')
title('conservation errors')
xlabel t
