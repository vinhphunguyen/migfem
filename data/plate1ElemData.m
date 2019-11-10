% rectangular plate in tension

% controlPts = zeros(2,2,2);
% 
% controlPts(1,1,:) = [0 0];
% controlPts(1,2,:) = [1 0];
% controlPts(2,1,:) = [0 1];
% controlPts(2,2,:) = [1 1];

controlPts=[0 0; 1 0;0 1; 1 1];
%controlPts=[0 0; 0 1;1 0; 1 1];
% knot vectors

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

weights = ones(1,4)';

p = 1;
q = 1;

noPtsX = 2;
noPtsY = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% controlPts=[0 0; 0 0.5; 0 1;
%             0.5 0;0.5 0.5;0.5 1;
%             1 0; 1 0.5; 1 1];
% 
% 
% % knot vectors
% 
% uKnot = [0 0 0 1 1 1];
% vKnot = [0 0 0 1 1 1];
% 
% weights = ones(1,9)';
% 
% p = 2;
% q = 2;
% 
% noPtsX = 3;
% noPtsY = 3;
