function [x,w]=gauleg(x1,x2,n)
%function computes Gaussian points
%inputs - Integral limits, x1 and x2
%         Number of gaussian points - n
%Outputs - x and W
%---------------------------------------------
eps=1e-15;
m=(n+1)/2;
xm=(1/2.)*(x2+x1);
xl=(1/2.)*(x2-x1);
pi=acos(-1.);

z1=0;

for i=1:m
    z=cos(pi*(i-0.25)/(n+0.5));
    p1=1.;
    p2=0.;
    
    while (abs(z-z1) > eps)
        p1=1.;
        p2=0.;
        for j=1:n
            p3=p2;
            p2=p1;
            p1=((2.*j-1.)*z*p2-(j-1.)*p3)/j;
        end
        pp=n*(z*p1-p2)/(z^2-1.);
        z1=z;
        z=z1-p1/pp;
    end
    x(i)=xm-x1*z;
    x(n+1-i)=xm+x1*z;
    w(i)=2.*x1/((1.-z*z)*pp*pp);
    w(n+1-i)=w(i);
end

x = x';
w = w' ;

end