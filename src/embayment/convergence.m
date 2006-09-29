% After updates 2006-07-27 with implicit draging
a = [25000 3.913296e-01 7.764259e-01;
     12500 1.569619e-01 2.372217e-01;
     6250 6.899825e-02 8.062357e-02;
     4166.67 4.392515e-02 4.530655e-02;
     2083.33 2.093577e-02 1.843607e-02;
     1388.89 1.373099e-02 1.135594e-02]
delta = a(:,1) * 1e-3;
u_error = a(:,2);
v_error = a(:,3);

clf
hold off
loglog(delta, u_error, 'rx--')
hold on
loglog(delta, v_error, 'g+-')
title('Refinement path for Macayeal Equations')
xlabel('\Delta x (km)')
ylabel('|u - u_{exact}| (m/a)')
legend('u', 'v')
