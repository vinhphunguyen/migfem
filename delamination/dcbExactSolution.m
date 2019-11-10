w  = 20;
h  = 1.5;
a0 = 30;

E   = 2.1e5;
G1c = 0.28;

I  = w*h^3/12;

% elastic branch
fileName = '~/code/jive/bezier/delamination/dcb3d/dcbExact1.dat';
u1 = linspace(0,5,20);
p1 = 1.5*E*I/a0^3*u1;

data = [u1' p1'];
save(fileName, 'data', '-ASCII');


% delamination branch
fileName = '~/code/jive/bezier/delamination/dcb3d/dcbExact2.dat';
u2 = linspace(1,10,20);
p2 = (w*G1c)^0.75*(E*I)^0.25*sqrt(2/3)*u2.^(-0.5);

data = [u2' p2'];
save(fileName, 'data', '-ASCII');

