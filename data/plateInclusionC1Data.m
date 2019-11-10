% the following is for the plate with a circular
% inclusion of unit radius.
% To model the material interface, knot with a multiplicity of
% 2 (=p) is used.
% Vinh Phu Nguyen
% Johns Hopkins University

tanPi8 = tan(pi/8);

controlPts        = [0 0;0 0;0 0;0 0; 
                     -0.3 0; -0.3 0.3*tanPi8;-0.3*tanPi8 0.3; 0 0.3;                     
                     -0.6 0;-0.6 tanPi8;-tanPi8 0.6; 0 0.6;
                     -1 0;-1 tanPi8;-tanPi8 1; 0 1;
                     -1.5 0;-1.5 1.5*tanPi8;-1.5*tanPi8 1.5; 0 1.5;
                     -3.1 0;-3.1 3.1*tanPi8;-3.1*tanPi8 3.1; 0 3.1;
                     -4 0;-4 4;-4 4;0 4];
                 
uKnot = [0 0 0 0.5 1 1 1];
vKnot = [0 0 0 0.15 0.3 0.3 0.7 1 1 1];

noPtsX = 4;
noPtsY = 7;

p     = 2;
q     = 2;

cont = 0.5*(1+1/sqrt(2));

weights = ones(noPtsX*noPtsY,1);

weights([6,7,10,11,14,15,18,19,22,23]) = cont;

    
