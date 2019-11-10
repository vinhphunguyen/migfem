function [stress disp]=exact_solution_timoshenkobeam(strPoint,E0,nu0,I,L,D,P)

x = strPoint(1,1);
y = strPoint(1,2);

stress(1) = (1/I)*P*(L-x)*y;
stress(2) = 0;
stress(3) = -0.5*(P/I)*(0.25*D^2 - y^2);

disp(1) = P*y/(6*E0*I)*((6*L-3*x)*x+(2+nu0)*(y^2-D^2/4));
disp(2) = -P/(6*E0*I)*(3*nu0*y^2*(L-x)+(4+5*nu0)*D^2*x/4+...
          (3*L-x)*x^2);