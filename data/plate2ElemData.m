% rectangular plate in tension

controlPts = [0 0; 0.3 0; 0.6 0; 1 0;
            %  0 0.3; 0.3 0.3; 0.6 0.3; 1 0.3;
            %  0 0.6; 0.3 0.6; 0.6 0.6; 1 0.6;
              0 1; 0.3 1; 0.6 1; 1 1];

% the following give negative Jacobian!!!
% controlPts = [0 0; 0 0.3;0 0.6; 0 1;
%               0.3 0;0.3 0.3; 0.3 0.6; 0.3 1;
%               0.6 0; 0.6 0.3;0.6 0.6;0.6 1;
%               1 0; 1 0.3; 1 0.6 ; 1 1];
          
% knot vectors

uKnot = [0 0 0 0.5 1 1 1];
vKnot = [0 0 1 1];

p = 2;
q = 1;

noPtsX = 4;
noPtsY = 2;

weights = ones(1,noPtsX*noPtsY)';